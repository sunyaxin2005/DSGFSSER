function [ W ] = Get_01Weight( X, L )
%GET_01WEIGHT Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(X);
K = zeros(N, N);
for ii = 1 : N
    nl = length(L(L == L(ii)));
    for jj = 1 : N
        if L(ii) == L(jj)
            K(ii, jj) = 1/nl;
        end
    end
end
W = K;
end
