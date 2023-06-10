function [ output_args ] = SPFS_LAR( X, Y)
%SPFS_LAR Summary of this function goes here
%   Detailed explanation goes here
k = size(Y, 2);
[N, D] = size(X);
XL = zeros(N * k, D * K);
for ii = 1 : D
    for jj = 1 : k
        XL((ii -1) * k + jj, N * (jj - 1) + 1 : N * (jj - 1) + k) = X(:, ii);
    end
end

end

