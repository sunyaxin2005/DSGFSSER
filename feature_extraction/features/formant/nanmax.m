function [f_max] = nanmax(data);
%
%   [f_max] = nanmax(data);
%
%Function which calculates the max (not NaN) of data containing
%NaN's.  NaN's are excluded completely from calculation.


[m,n] = size(data);

for index = 1:n;
not_nans = find(isnan(data(:,index)) == 0);
    if length(not_nans) > 0;
        f_max(index) = max(data(not_nans,index));
    else
        f_max(index) = NaN;
    end
end