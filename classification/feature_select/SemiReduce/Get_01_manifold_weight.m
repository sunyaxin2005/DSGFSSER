function [ W ] = Get_01_manifold_weight(  X, L )
%GET_01_MANIFOLD_WEIGHT Summary of this function goes here
%   Detailed explanation goes here
k = 7;
[N, D] = size(X);
K = zeros(N, N);
distance = pdist2(X, X);
[sorted,index] = sort(distance);
[ W ] = get_heat_kernel( X );
for ii = 1 : N
    for jj = 1 : k
        K(ii, index(ii, jj)) = K(ii, index(ii, jj)) + W(ii, index(ii, jj));
    end
end
for ii = 1 : N
    nl = length(L(L == L(ii)));
    for jj = 1 : N
        if L(ii) == L(jj) && L(ii) ~= -1
            K(ii, jj) = 1;
        end
    end
    K(ii, :) = K(ii, :)./sum(K(ii, :));
end

W = K + K';
end

