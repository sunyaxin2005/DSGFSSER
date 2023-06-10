function [ W ] = Get_MI_Weight( X)
%GET_MI_WEIGHT Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(X);
W = zeros(D, D);
for ii = 1 : D
   for jj = 1 : D
        W(ii, jj) = mi_right(X(:, ii), X(:, jj), 32);
    end
end
Z = W;
Z(:, :) = 0;
for ii = 1 : N
    inx = neighborhood(:, ii);
    Z(:, ii) = W(:, ii)./sum(W(:, ii)); 
end
%W = max(Z, Z');
W = Z + Z' - Z * Z';
end

