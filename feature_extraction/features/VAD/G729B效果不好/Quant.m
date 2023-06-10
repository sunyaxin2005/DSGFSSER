function Index = Quant (x, Xq, QType)
% Quantize a vector of samples
%
% x - vector of input samples
% Xq - Quantizer decision levels (Nreg-1 levels, where Nref is the number
%    of quantizer intervals).
% QType - Type of quantizer (optional, default 1). These quantizer types
%    differ at the end points of the intervals. For Type 1, if an input
%    value is equal to the upper end of an interval, it is considered to
%    lie in the interval. For Type 2, if an input value is equal to the
%    lower end of an interval, it is consider to lie in the interval.
%
% Index - Output index (0 <= Index < Nreg)
%
%  This function returns the index of the quantizer region corresponding to
%  a given input value. The quantizer is specified by an array of quantizer
%  decision levels. A binary search is used to determine the quantizer region.
%
%  The index value takes on values from 0 to Nreg-1, where Nreg is the number
%  of quantizer regions. If number of regions is equal to one (Xq is empty),
%  the index is set to zero. Otherwise, the index is determined as shown in
%  the following table. The number of decision levels is one less than the
%  number of regions.
%    Index          Interval QType 1              Interval QType 2
%     0             -Inf < x <= Xq(1)            -Inf <= x < Xq(1)
%     1            Xq(1) < x <= Xq(2)           Xq(1) <= x < Xq(2)
%     k            Xq(k) < x <= Xq(k+1)         Xq(k) <= x < Xq(k+1)
%   Nlev-1    Xq(Nlev-1) < x <= Inf        Xq(Nlev-1) <= x < Inf

% For Nlev = 256, and random uniform data, this routine is about 2 times
% faster than MATLAB's quantiz

Nreg = length(Xq) + 1;
Nval = length(x);
Index = zeros(size(x));

if (nargin <= 2)
  QType = 0;
end

if (QType == 2)

% ===== =====
  for (n = 1:Nval)
    iL = 0;
    iU = Nreg;

% Isolate the interval
    while (iU > iL + 1)
      i = fix((iL + iU) / 2);
      if (x(n) < Xq(i))
        iU = i;
      else
        iL = i;
      end
    end

    Index(n) = iL;
  end

else

% ===== =====
  for (n = 1:Nval)
    iL = 0;
    iU = Nreg;

% Isolate the interval
    while (iU > iL + 1)
      i = fix((iL + iU) / 2);
      if (x(n) <= Xq(i))
        iU = i;
      else
        iL = i;
      end
    end

    Index(n) = iL;
  end
  
% ====== =====
  
end

return
