function [f_min] = nanmin(data);
%
%   [f_min] = nanmin(data);
%
%Function which calculates the min (not NaN) of data containing
%NaN's.  NaN's are excluded completely from calculation.


[m,n] = size(data);

for index = 1:n;
not_nans = find(isnan(data(:,index)) == 0);
    if length(not_nans) > 0;
        f_min(index) = min(data(not_nans,index));
    else
        f_min(index) = NaN;
    end
end