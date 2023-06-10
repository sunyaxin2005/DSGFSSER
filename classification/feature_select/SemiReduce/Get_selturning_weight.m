function [ W ] = Get_selturning_weight( X, K )
%GET_SELTURNING_WEIGHT Summary of this function goes here
%   Detailed explanation goes here
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
W = exp(-(distance.^2)./sita);

end

