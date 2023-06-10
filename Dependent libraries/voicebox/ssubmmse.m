function [ss,gg,tt,ff,zo]=ssubmmse(si,fsz,pp)
%SSUBMMSE performs speech enhancement using mmse estimate of spectral amplitude or log amplitude [SS,ZO]=(S,FSZ,P)
%
% Usage: y=ssubmmse(x,fs);   % enhance the speech using default parameters
%
% Inputs:
%   si      input speech signal
%   fsz     sample frequency in Hz
%           Alternatively, the input state from a previous call (see below)
%   pp      algorithm parameters [optional]
%
% Outputs:
%   ss        output enhanced speech
%   gg(t,f,i) selected time-frequency values (see pp.tf below)
%   tt        centre of frames (in seconds)
%   ff        centre of frequency bins (in Hz)
%   zo        output state (or the 2nd argument if gg,tt,ff are omitted)
%
% The algorithm operation is controlled by a small number of parameters:
%
%        pp.of          % overlap factor = (fft length)/(frame increment) [2]
%        pp.ti          % desired frame increment [0.016 seconds]
%        pp.ri          % set to 1 to round ti to the nearest power of 2 samples [0]
%        pp.ta          % time const for smoothing SNR estimate [0.396 seconds]
%        pp.gx          % maximum posterior SNR as a power ratio [1000 = +30dB]
%        pp.gn          % min posterior SNR as a power ratio when estimating prior SNR [1 = 0dB]
%        pp.gz          % min posterior SNR as a power ratio [0.001 = -30dB]
%        pp.xn          % minimum prior SNR [0]
%        pp.xb          % bias compensation factor for prior SNR [1]
%        pp.lg          % MMSE target: 0=amplitude, 1=log amplitude [1]
%        pp.ne          % noise estimation: 0=min statistics, 1=MMSE [0]
%        pp.bt          % threshold for binary gain or -1 for continuous gain [-1]
%        pp.mx          % input mixture gain [0]
%        pp.rf          % round output signal to an exact number of frames [0]
%        pp.tf          % selects time-frequency planes to output in the gg() variable ['g']
%                           'i' = input power spectrum
%                           'I' = input complex spectrum
%                           'n' = noise power spectrum
%                           'z' = "posterior" SNR (i.e. (S+N)/N )
%                           'x' = "prior" SNR (i.e. S/N )
%                           'g' = gain
%                           'o' = output power spectrum
%                           'O' = output complex spectrum
%
% The applied gain is mx+(1-mx)*optgain where optgain is calculated according to [1] or [2].
% If pp.bt>=0 then optgain is first thresholded with pp.bt to produce a binary gain 0 or 1.
%
% The default parameters implement the original algorithm in [1,2].
%
% Several parameters relate to the estimation of xi, the so-called "prior SNR",
%
%             xi=max(a*pp.xb*xu+(1-a)*max(gami-1,pp.gn-1),pp.xn);
%
% This is estimated as a smoothed version of gami, the "posterior SNR"
% which is the noisy speech power divided by the noise power. This is
% clipped to a min of (pp.gn-1), smoothed using a factor "a" which corresponds to a
% time-constant of pp.ta and then clipped to a minimum of pp.xn. The
% previous value is taken to be pp.xb*xu where xu is the ratio of the
% estimated speech amplitude squared to the noise power.
%
% In addition it is possible to specify parameters for the noise estimation algorithm
% which implements reference [3] or [7] according to the setting of pp.ne
% 
% Minimum statistics noise estimate [3]: pp.ne=0 
%        pp.taca      % (11): smoothing time constant for alpha_c [0.0449 seconds]
%        pp.tamax     % (3): max smoothing time constant [0.392 seconds]
%        pp.taminh    % (3): min smoothing time constant (upper limit) [0.0133 seconds]
%        pp.tpfall    % (12): time constant for P to fall [0.064 seconds]
%        pp.tbmax     % (20): max smoothing time constant [0.0717 seconds]
%        pp.qeqmin    % (23): minimum value of Qeq [2]
%        pp.qeqmax    % max value of Qeq per frame [14]
%        pp.av        % (23)+13 lines: fudge factor for bc calculation  [2.12]
%        pp.td        % time to take minimum over [1.536 seconds]
%        pp.nu        % number of subwindows to use [3]
%        pp.qith      % Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf ]
%        pp.nsmdb     % corresponding noise slope thresholds in dB/second   [47 31.4 15.7 4.1]
%
% MMSE noise estimate [7]: pp.ne=1 
%        pp.tax      % smoothing time constant for noise power estimate [0.0717 seconds](8)
%        pp.tap      % smoothing time constant for smoothed speech prob [0.152 seconds](23)
%        pp.psthr    % threshold for smoothed speech probability [0.99] (24)
%        pp.pnsaf    % noise probability safety value [0.01] (24)
%        pp.pspri    % prior speech probability [0.5] (18)
%        pp.asnr     % active SNR in dB [15] (18)
%        pp.psini    % initial speech probability [0.5] (23)
%        pp.tavini   % assumed speech absent time at start [0.064 seconds]
%
% If convenient, you can call specsub in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) y=ssubmmse(s,fs);
%
%                   (b) [y1,z]=ssubmmse(s(1:1000),fs);
%                       [y2,z]=ssubmmse(s(1001:2000),z);
%                       y3=ssubmmse(s(2001:end),z);
%                       y=[y1; y2; y3];
%
% If the number of output arguments is either 2 or 5, the last partial frame of samples will
% be retained for overlap adding with the output from the next call to ssubmmse().
%
% See also specsub() for an alternative gain function
%
% Refs:
%    [1] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
%    [2] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum mean-square error log-spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 33(2):443-445, Apr 1985
%    [3] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
%    [4] O. Cappe.
%        Elimination of the musical noise phenomenon with the ephraim and malah noise suppressor.
%        IEEE Trans Speech Audio Processing, 2 (2): 345�349, Apr. 1994. doi: 10.1109/89.279283.
%    [5] J. Erkelens, J. Jensen, and R. Heusdens.
%        A data-driven approach to optimizing spectral speech enhancement methods for various error criteria.
%        Speech Communication, 49: 530�541, 2007. doi: 10.1016/j.specom.2006.06.012.
%    [6] R. Martin.
%        Statistical methods for the enhancement of noisy speech.
%        In J. Benesty, S. Makino, and J. Chen, editors,
%        Speech Enhancement, chapter 3, pages 43�64. Springer-Verlag, 2005.
%    [7] Gerkmann, T. & Hendriks, R. C.
%        Unbiased MMSE-Based Noise Power Estimation With Low Complexity and Low Tracking Delay
%        IEEE Trans Audio, Speech, Language Processing, 2012, 20, 1383-1393

% Bugs/suggestions:
%   (1) sort out behaviour when si() is a matrix rather than a vector
%
%      Copyright (C) Mike Brookes 2004-2011
%      Version: $Id: ssubmmse.m 1720 2012-03-31 17:17:31Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(si)>length(si)
    error('Input speech signal must be a vector not a matrix');
end
if isstruct(fsz)
    fs=fsz.fs;
    qq=fsz.qq;
    qp=fsz.qp;
    ze=fsz.ze;
    s=zeros(length(fsz.si)+length(si(:)),1); % allocate space for speech
    s(1:length(fsz.si))=fsz.si;
    s(length(fsz.si)+1:end)=si(:);
else
    fs=fsz;     % sample frequency
    s=si(:);
    % default algorithm constants

    qq.of=2;        % overlap factor = (fft length)/(frame increment)
    qq.ti=16e-3;    % desired frame increment (16 ms)
    qq.ri=0;        % round ni to the nearest power of 2
    qq.ta=0.396;    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [1]
    qq.gx=1000;     % maximum posterior SNR = 30dB
    qq.gn=1;        % min posterior SNR as a power ratio when estimating prior SNR [1]
    qq.gz=0.001;    % min posterior SNR as a power ratio [0.001 = -30dB]
    qq.xn=0;        % minimum prior SNR = -Inf dB
    qq.xb=1;        % bias compensation factor for prior SNR [1]
    qq.lg=1;        % use log-domain estimator by default
    qq.ne=0;        % noise estimation: 0=min statistics, 1=MMSE [0]
    qq.bt=-1;       % suppress binary masking
    qq.mx=0;        % no input mixing
    qq.tf='g';      % output the gain time-frequency plane by default
    qq.rf=0;
    if nargin>=3 && ~isempty(pp)
        qp=pp;      % save for estnoisem call
        qqn=fieldnames(qq);
        for i=1:length(qqn)
            if isfield(pp,qqn{i})
                qq.(qqn{i})=pp.(qqn{i});
            end
        end
    else
        qp=struct;  % make an empty structure
    end
end
% derived algorithm constants
if qq.ri
    ni=pow2(nextpow2(ti*fs*sqrt(0.5)));
else
    ni=round(qq.ti*fs);    % frame increment in samples
end
tinc=ni/fs;         % true frame increment time
a=exp(-tinc/qq.ta); % SNR smoothing coefficient
gx=qq.gx;           % max posterior SNR as a power ratio
gz=qq.gz;           % min posterior SNR as a power ratio
kk=sqrt(2*pi);      % sqrt(8)*Gamma(1.5) - required constant
xn=qq.xn;           % floor for prior SNR, xi
ne=qq.ne;           % noise estimation: 0=min statistics, 1=MMSE [0]
gn1=max(qq.gn-1,0); % floor for posterior SNR when estimating prior SNR
xb=qq.xb;
tf=qq.tf;
rf=qq.rf || nargout==2 || nargout==5;            % round down to an exact number of frames

% calculate power spectrum in frames

no=round(qq.of);                  	% integer overlap factor
nf=ni*no;                           % fft length
w=sqrt(hamming(nf+1))'; w(end)=[];  % for now always use sqrt hamming window
w=w/sqrt(sum(w(1:ni:nf).^2));       % normalize to give overall gain of 1
if rf>0
    rfm='';                         % truncated input to an exact number of frames
else
    rfm='r';
end
[y,tt]=enframe(s,w,ni,rfm);
tt=tt/fs;                           % frame times
yf=rfft(y,nf,2);
yp=yf.*conj(yf);                    % power spectrum of input speech
[nr,nf2]=size(yp);                  % number of frames
ff=(0:nf2-1)*fs/nf;
if isstruct(fsz)
    if ne>0
        [dp,ze]=estnoiseg(yp,ze);       % estimate the noise using MMSE
    else
        [dp,ze]=estnoisem(yp,ze);       % estimate the noise using minimum statistics
    end
    ssv=fsz.ssv;
    xu=fsz.xu;                      % saved unsmoothed SNR
else
    if ne>0
        [dp,ze]=estnoiseg(yp,tinc,qp);	% estimate the noise using MMSE
    else
        [dp,ze]=estnoisem(yp,tinc,qp);	% estimate the noise using minimum statistics
    end
    ssv=zeros(ni*(no-1),1);        	% dummy saved overlap
    xu=1;                           % dummy unsmoothed SNR from previous frame
end
if ~nr                             	% no data frames
    ss=[];
    gg=[];
else
    gam=max(min(yp./dp,gx),gz);     % gamma = posterior SNR
    g=zeros(nr,nf2);                % create space for gain matrix
    x=zeros(nr,nf2);                % create space for prior SNR
    if qq.lg                        % use log domain estimator
        for i=1:nr
            gami=gam(i,:);
            xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
            xir=xi./(1+xi);
            gi=xir.*exp(0.5*expint(xir.*gami));
            g(i,:)=gi;             	% save gain for later
            x(i,:)=xi;              % save prior SNR
            xu=gami.*gi.^2;         % unsmoothed prior SNR
        end
    else
        for i=1:nr
            gami=gam(i,:);
            xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
            v=0.5*xi.*gami./(1+xi);	% note that this is 0.5*vk in [1]
            gi=(0.277+2*v)./gami; 	% accurate to 0.02 dB for v>0.5
            mv=v<0.5;
            if any(mv)
                vmv=v(mv);
                gi(mv)=kk*sqrt(vmv).*((0.5+vmv).*besseli(0,vmv)+vmv.*besseli(1,vmv))./(gami(mv).*exp(vmv));
            end
            g(i,:)=gi;              % save gain for later
            x(i,:)=xi;              % save prior SNR
            xu=gami.*gi.^2;         % unsmoothed prior SNR
        end
    end
    if qq.bt>=0
        g=g>qq.bt;
    end
    g=qq.mx+(1-qq.mx)*g;                    % mix in some of the input
    se=(irfft((yf.*g).',nf).').*repmat(w,nr,1);     % inverse dft and apply output window
    ss=zeros(ni*(nr+no-1),no);                      % space for overlapped output speech
    ss(1:ni*(no-1),end)=ssv;
    for i=1:no
        nm=nf*(1+floor((nr-i)/no));         % number of samples in this set
        ss(1+(i-1)*ni:nm+(i-1)*ni,i)=reshape(se(i:no:nr,:)',nm,1);
    end
    ss=sum(ss,2);
    if nargout>2 && ~isempty(tf)
        gg=zeros(nr,nf2,length(tf));  % make space
        for i=1:length(tf)
            switch tf(i)
                case 'i'            % 'i' = input power spectrum
                    gg(:,:,i)=yp;
                case 'I'            % 'i' = input power spectrum
                    gg(:,:,i)=yf;
                case 'n'            % 'n' = noise power spectrum
                    gg(:,:,i)=dp;
                case 'z'            % 'z' = posterior SNR (i.e. (S+N)/N )
                    gg(:,:,i)=gam;
                case 'x'            % 'x' = prior SNR
                    gg(:,:,i)=x;
                case 'g'            % 'g' = gain
                    gg(:,:,i)=g;
                case 'o'            % 'o' = output power spectrum
                    gg(:,:,i)=yp.*g.^2;
                case 'O'            % 'o' = output power spectrum
                    gg(:,:,i)=yf.*g;
            end
        end
    end
end
if nargout==2 || nargout==5
    if nr
        zo.ssv=ss(end-ni*(no-1)+1:end);    % save the output tail for next time
        ss(end-ni*(no-1)+1:end)=[];        % only output the frames that are completed
    else
        zo.ssv=ssv;  %
    end
    zo.si=s(length(ss)+1:end);      % save the tail end of the input speech signal
    zo.fs=fs;                       % save sample frequency
    zo.qq=qq;                       % save local parameters
    zo.qp=qp;                       % save estnoisem parameters
    zo.ze=ze;                       % save state of noise estimation
    zo.xu=xu;
    if nargout==2
        gg=zo;                      % 2nd of two arguments is zo
    end
elseif rf==0
    ss=ss(1:length(s));             % trim to the correct length if not an exact number of frames
end
if ~nargout && nr>0
    ffax=ff/1000;
    ax=zeros(4,1);
    ax(1)=subplot(223);
    imagesc(tt,ffax,20*log10(g)');
    colorbar;
    axis('xy');
    title(sprintf('Filter Gain (dB): ta=%.2g',qq.ta));
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(2)=subplot(222);
    imagesc(tt,ffax,10*log10(yp)');
    colorbar;
    axis('xy');
    title('Noisy Speech (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(3)=subplot(224);
    imagesc(tt,ffax,10*log10(yp.*g.^2)');
    colorbar;
    axis('xy');
    title('Enhanced Speech (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');

    ax(4)=subplot(221);
    imagesc(tt,ffax,10*log10(dp)');
    colorbar;
    axis('xy');
    title('Noise Estimate (dB)');
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    linkaxes(ax);
end