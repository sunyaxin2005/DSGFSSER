function [ window_res ] = hanning_window( x )
%HANNING_WINDOW Summary of this function goes here
%   Detailed explanation goes here
[rows cols] = size(x);
h = hamming(rows);
window_res = zeros(rows, cols);
for n = 1 : cols
    for m = 1 : rows
        window_res(m, n) = h(m) * x(m, n);
    end
end

end

