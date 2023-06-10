function DTXPar = InitDTXPar()

% Constants
DTXPar.M = 10;
DTXPar.NFrmAvgAcf = 3;
DTXPar.NFrmCurAcf = 2;
DTXPar.NSIDMin = 3;
DTXPar.Thresh1 = 1.1481628;
DTXPar.Thresh2 = 1.0966466;
DTXPar.ThreshEndB = 2;
DTXPar.NEn = 2;

DTXPar.N = 240;    % window size
DTXPar.LA = 40;    % Look-ahead
DTXPar.NF = 80;    % Frame size

LWmem = DTXPar.N - DTXPar.NF;
DTXPar.Wmem = zeros(LWmem, 1);

% HPFilt is a HPF that is used to preprocess the signal applied to the DTX module.
DTXPar.HPFilt.b = [ 0.92727435, -1.8544941,  0.92727435 ];
DTXPar.HPFilt.a = [ 1,          -1.9059465,  0.91140240 ];
DTXPar.HPFilt.Mem = [];

LA = DTXPar.LA;
LB = DTXPar.N - DTXPar.LA;
DTXPar.Window = [0.54 - 0.46*cos(2*pi*(0:LB-1)'/(2*LB-1));
                 cos(2*pi*(0:LA-1)'/(4*LA-1))];

% The energy scaling factor G.729C+ is documented as being composed of
% several terms
%   fact[n] = fact_ener / (n x NF x nbAcf)
%   fact = [0.00078125, 0.000390625]  for n = 1 and 2.
% By taking the mean of the energies rather than the sum, we eliminate the
% n term. By taking the mean rather than the sum of correlations, we
% eliminate the nbAcf term. The remaining term is a fudge factor which
% includes the window scaling. From these considerations, fact_ener = 1/8.
fact_ener = 1/8;
DTXPar.EnQ.Scale = fact_ener / DTXPar.NF;

% Energy quantizer (dB) units
Xq = [   -8   -2   2   6  10  14  17  19  21  23  25  27  29  31  33 ...
      35  37  39  41  43  45  47  49  51  53  55  57  59  61  63  65 ];
Yq = [-12  -4   0   4   8  12  16  18  20  22  24  26  28  30  32  34 ...
       36  38  40  42  44  46  48  50  52  54  56  58  60  62  64  66];
DTXPar.EnQ.Xq = Xq;
DTXPar.EnQ.Yq = Yq;

DTXPar.VadPrev = 1;    % Previous VAD value
DTXPar.FrmSID = 0;     % Number of frames since the last SID frame
DTXPar.FrmCurAcf = 0;  % Counter to signal update to Acf (during active)
DTXPar.Acf = zeros(DTXPar.M+1, DTXPar.NFrmCurAcf);
DTXPar.AvgAcf = zeros(DTXPar.M+1, DTXPar.NFrmAvgAcf);
DTXPar.En = zeros(1, DTXPar.NEn);
DTXPar.nE = 0;         % Number of energies available
DTXPar.EndBPrev = -Inf;
DTXPar.EndB = -Inf;
DTXPar.ChangeFlag = 0; % Signal a change in noise

return
