function [ mfcc_res ] = wsf( x, fs, frame_long, inc, alg_prog )
%WSF Summary of this function goes here
%   Detailed explanation goes here
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[ sample_num frame_num] = size(f);
[ b ] = mel_fft( x, fs, frame_long, inc, alg_prog);
for ii = 1 : 1:frame_num
    ar = lpcauto(f(:, ii), 10);
    [row, p] = size(ar);
    hlpc2lsp = dsp.LPCToLSF;
    y = step(hlpc2lsp, ar');
    p = p - 1; 
    w = zeros(1, p);
    w(1) = 1/(y(2) - y(1));
    for jj = 2 : p-1
        w(jj) = 1/(y(jj) - y(jj - 1)) + 1/(y(jj + 1) - y(jj));
    end
    w(p) = 1/(y(p) - y(p-1));
    w = w./sum(w);
    w = repmat(w, 4, 1);
    w = w(:)';
    gauss = fspecial('gaussian', [1 4], 1);
    w = imfilter(w, gauss);
    b(ii, :) = b(ii, :).*w;
end
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: 13)];
end
end

function [ b ] = mel_fft( s, fs, frame_long, inc, alg_prog)
%MEL_FFT Summary of this function goes here
%   Detailed explanation goes here
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
%b = log(b);
end