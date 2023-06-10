function [ b ] = mel_fft( s, fs, frame_long, inc, alg_prog)
%MEL_FFT Summary of this function goes here
%   Detailed explanation goes here
%x,   speech data
%fs, the sampling rate
%frame_long, Hamming window width
%inc, Hamming window shift
%alg_prog, the number of Mel filters  
[ns1,ns2]=size(s);
winlen = frame_long;
ninc = inc;
tinc=[0 -Inf Inf];
fs(2)=1.0/fs(1);
nfrq = alg_prog;             
win=0.54+0.46*cos((1-winlen:2:winlen)*pi/winlen);
fftlen=pow2(nextpow2(4*winlen));        % enough oversampling to get good interpolation
win=win/sqrt(sum(win.^2));              % ensure window squared sums to unity
ix1=max(round((tinc(2)-fs(2))*fs(1)-(winlen-3)/2),1); % first sample required
ix2=min(ceil((tinc(3)-fs(2))*fs(1)+(winlen+1)/2),ns1); % last sample required
[sf,t]=enframe(s(ix1:ix2),win,ninc);
t=fs(2)+(t+ix1-2)/fs(1);                         % time axis
b=rfft(sf,fftlen,2);
b=b.*conj(b)*2/fs(1);          % Power per Hz
b(:,1)=b(:,1)*0.5;   % correct for no negative zero frequency to double the power
b(:,end)=b(:,end)*0.5;   % correct for no negative nyquist frequency to double the power
fb=(0:fftlen/2)*fs(1)/fftlen; % fft bin frequencies
fftfs=fs(1);
nfr=numel(t);  
b=b.*repmat((700+fb)*log(1+1000/700)/1000,nfr,1); 
b=b*filtbankm(nfrq,fftlen,fftfs,0,4000,'m')';
b = log(b);
end

