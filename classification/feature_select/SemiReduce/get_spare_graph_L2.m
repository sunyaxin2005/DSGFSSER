function [ S ] = get_spare_graph_L2( X)
%SPARE_GRAPH Summary of this function goes here
%   Detailed explanation goes here
 %       X = pdist2(X, X);
%         [sorted,index] = sort(distance);
%         sita = sorted(K+1, :);
%         sita = sita' *  sita;
%         X = exp(-(distance.^2)./sita);
n = size(X, 1);
m = size(X, 2);
X = X';

labda = 0.1;
I = eye(n,n);
Q = (X'*X + labda * I);

e = zeros(n ,1);
S = zeros(n, n);
for ii = 1 : n
    zi = X' * X(:, ii);
    ei = e;
    ei(ii) = 1;
    si = Q * (zi - (ei' * Q * zi * ei) /(ei' * Q * ei)); 
    S(ii, :) = si;
end


end
