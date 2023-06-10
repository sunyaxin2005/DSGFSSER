function [ vad_res ] = VAD_main( x, fs )
%VAD_MAIN Summary of this function goes here
%   Detailed explanation goes here

Ns = length(x);
NsFrame = 80;       % Frame size (10ms)
NsLA = 40;          % Look-ahead size for VAD / DTX
NsWin = 240;
NFrame = floor(Ns/NsFrame);

x = [x; zeros(NsFrame, 1)];
VADPar = InitVADPar(fs);
vad_res = [];
for (k = 0:NFrame-1)

  ist = k*NsFrame + 1;
  ifn = ist + NsFrame - 1;       % New data limits

  % VAD / DTX
  x_new = x(ist:ifn);
  [Ivd, VADPar] = VAD(x_new, VADPar);
  vad_res = [vad_res; Ivd];
end

end

