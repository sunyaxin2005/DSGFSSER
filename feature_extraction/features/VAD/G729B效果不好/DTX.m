function [Ftype, DTXPar] = DTX (x_new, VAD, DTXPar)
% DTX module for G.729B / G.729C+

% Filter new data (HP filter)
[x_new_hp, DTXPar.HPFilt.Mem] = filter(DTXPar.HPFilt.b, DTXPar.HPFilt.a, ...
                                       32768 * x_new, DTXPar.HPFilt.Mem);

% Append new filtered data to filter memory
xwin = [DTXPar.Wmem; x_new_hp];

% Compute autocorrelation
r = acorr(DTXPar.Window .* xwin, DTXPar.M + 1);
DTXPar.Acf = [r, DTXPar.Acf(:,1:end-1)];    % Push new correlation values

% Notes:
% In G.729 DTX, the correlations are summed, not averaged. Later on during
% the quantization procedure, a normalization of the energy compensates
% for several factors including frame length and the number of correlation
% vectors which are used in the averages. Here we use averages.
%
% In G.729 DTX, the energies are summed and then normalized to take out
% the effects of summing.

% P. Kabal 2008-04-03

% ----- -----
% Frame type
if (VAD == 1)  
  % Active frame
  Ftype = 1;

else
  % Inactive frame

% New filter for this frame
  [CurFilt, DTXPar.En] = GenCurFilt(DTXPar);

% Determine if an SID frame needs to be transmitted
  [Ftype, DTXPar] = CheckSID (CurFilt, DTXPar);

end

if (Ftype == 2)
  % Update parameters as per an SID frame
  DTXPar = ProcSID(CurFilt, DTXPar);
end

DTXPar.VADPrev = VAD;

% Periodically update the short-term mean correlations
DTXPar.FrmCurAcf = DTXPar.FrmCurAcf + 1;
if (DTXPar.FrmCurAcf == DTXPar.NFrmCurAcf)

% Update stack of mean correlation vectors
  DTXPar.FrmCurAcf = 0;          % Reset counter
  SAcf = mean(DTXPar.Acf, 2);    % Mean along rows
  DTXPar.AvgAcf = [SAcf, DTXPar.AvgAcf(:,1:end-1)];

end

% Update the analysis window memory
ist = DTXPar.NF + 1;
DTXPar.Wmem = xwin(ist:end);

return

% ======== ========
function [CurFilt, En] = GenCurFilt (DTXPar)
% Generate a new filter based on the short-term mean
% and update the vector of frame energies

% Current correlation is the short-term mean
MAcf = mean(DTXPar.Acf, 2);

% Current filter
CurFilt = Acf2Filt(MAcf);

% Update the list of energies
En = [CurFilt.En, DTXPar.En(1:end-1)];

return

% ======== ========
function [Ftype, DTXPar] = CheckSID (CurFilt, DTXPar)
% Process an inactive frame
% - Check that the current filter and energy of this noise frame has not
%   changed significantly from the filter and energy in the SID frame

if (DTXPar.VADPrev == 1)
  % First SID frame
  Ftype = 2;
  DTXPar.FrmSID = 0;      % Reset SID counter
  DTXPar.nE = 1;
  DTXPar.EndB = LogMeanEn(DTXPar.En(1:DTXPar.nE), DTXPar.EnQ);

else
  % Continuing inactive frame

  DTXPar.nE = min(DTXPar.nE + 1, DTXPar.NEn);
  DTXPar.EndB = LogMeanEn(DTXPar.En(1:DTXPar.nE), DTXPar.EnQ);

% Check for a change in the LP model
  if (distIS(CurFilt.Acf, DTXPar.SIDFilt.Acf) >= DTXPar.Thresh1)
    DTXPar.ChangeFlag = 1;
  end

% In G.729C+, it is the quantized energies that are used to check for
% a change in noise energy. These quantized energies are on a dB scale).
  if (abs(DTXPar.EnSIDdB - DTXPar.EndB) > DTXPar.ThreshEndB)
    DTXPar.ChangeFlag = DTXPar.ChangeFlag + 2;
  end

% Check if we should classify this inactive frame as an SID frame
  DTXPar.FrmSID = DTXPar.FrmSID + 1;
  if (DTXPar.FrmSID < DTXPar.NSIDMin)
    Ftype = 0;          % SID frames cannot be too close together
  else
    if (DTXPar.ChangeFlag ~= 0)
      Ftype = 2;        % Change in noise detected
    else
      Ftype = 0;
    end
    DTXPar.FrmSID = DTXPar.NSIDMin;  % Clamp at NSIDMin
  end

end

return

% ======== ========
function DTXPar = ProcSID (CurFilt, DTXPar)

DTXPar.FrmSID = 0;      % Reset the SID counter
DTXPar.ChangeFlag = 0;

% Find the long term average filter
AvgFilt = Acf2Filt(mean(DTXPar.AvgAcf, 2));

% Use the average filter if it is close enough
if (distIS(CurFilt.Acf, AvgFilt.Acf) < DTXPar.Thresh2)
  DTXPar.SIDFilt = AvgFilt;
else
  DTXPar.SIDFilt = CurFilt;
end

DTXPar.EnSIDdB = DTXPar.EndB;

return

% ======== ========
function Filt = Acf2Filt (Acf)

Filt.Acf = Acf;
[Filt.A, Filt.En] = ac2poly(Acf);

return

% ======== ========
function EndBQ = LogMeanEn (En, Quant)

SEn = mean(En) * Quant.Scale;
EndB = 10 * log10(SEn);
index = Quantz(EndB, Quant.Xq);
EndBQ = Quant.Yq(index+1);

return

% ======== ========
% Compute Itakura-Saito distance
function d = distIS (AcfRef, Acf)

global FIDx

% The Itakura-Saito distance measure is
%   d = (A'*RR*A) / (AR'*RR*AR),
% where RR is a Toeplitz matrix formed from the correlation vector.
% The denominator is the energy returned by the Levinson recursion. It
% represents the mean-square error for a predictor AR derived from the
% correlation values RR. The numerator represents the mean-square error
% for a (non-necessarily optimal) predictor acting on a signal with
% correlation RR.
%
% The numerator can be recast as a convolution. Let
%   RAA[k] = SUM_i A[i] A[i+k] = RAA[-k]
% Then
%   Ryy[n] = SUM_m RxxR[m] RAA[m+n].
% The numerator of the IS distance is Ryy[0].
%
% In G.729 DTX, the numerator is computed in this way. The numerator is
% evaluated using RAA (an input) and RxxR (another input). The denominator
% is directly an input (derived from the Levinson recursion). The
% computation of the convolution is reduced by taking into account the
% symmetry of the correlation values.
%
% The G.729 DTX routine returns a flag if the IS distortion exceeds a
% a threshold. Note that d is always larger than one. The comparison used
% in G.729 DTX is
%   A'*RR*A : AR'*RR*AR * Threshold

A = ac2poly(Acf);
A = A(:);
[ARef, ERef] = ac2poly(AcfRef);

R = toeplitz(AcfRef);

d = (A' * R * A) / ERef;

return

% -----------------------------
function rxx = acorr (x, Nt)

Nx = length (x);
N = Nt;
if (Nt > Nx)
  N = Nx;
end

rxx = zeros(Nt, 1);
for (i = 0:N-1)
  Nv = Nx - i;
  rxx(i+1) = x(1:Nv)' * x(i+1:i+Nv);
end

return

% ----- -----
function Index = Quantz (x, Xq)

Nreg = length(Xq) + 1;
Index = zeros(size(x));

[Temp, Index] = histc(x, Xq);
IU = (x > Xq(end));
Index(IU) = Nreg-1;

return