function [ W ] = GetSPPDW( X, k )
%GETSPPDW Summary of this function goes here
%   Detailed explanation goes here
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(k, :);
sita = sita' *  sita;
W = exp(-(distance.^2)./1);

end

