function [ W ] = get_heat_kernel( X )
%GET_HEAT_KERNEL Summary of this function goes here
%   Detailed explanation goes here
distance = pdist2(X, X);
sita = mean(mean(distance));
W = exp(-(distance.^2)./(sita * sita));
diagA = diag(ones(size(W, 2), 1));
W = W - diagA;
end

