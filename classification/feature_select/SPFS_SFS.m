function [ cLbest ] = SPFS_SFS( X, L, sel_dim )
%SPFS_SFS Summary of this function goes here
%   Detailed explanation goes here
%X: N * D
%L: 类别标签
%sel_dim:待选择的维度
if isempty(L)
    K = compute_K_ulsur(X);
else
    K = compute_K_sur(X, L);
end
R = K;
cLbest = [];
[N, D] = size(X);
pre_compute_values = zeros(D, N, N);
for ii = 1 : D
    pre_compute_values(ii) = X(:, ii)' *  X(:, ii);
end

for ii = 1 : sel_dim
    [res] = col_rese(pre_compute_values, R, cLbest);
    [C, I] = min(res);
    temp = zeros(N, N);
    temp(:, :) = pre_compute_values(ii, :, :);
    if norm(R - temp, 2) > norm(R, 2)
        break;
    else
        cLbest = [cLbest;I(1)];
        R = R - X(:, I(1))' *  X(:, I(1));
    end
end
end

function [res] = col_rese(pre_compute_values, R, cLbest)
[D, N, ~] = size(pre_compute_values);
res = zeros (1, D);
for ii = 1 : D
    indices = find(cLbest, ii, 'first');
    if isempty(indices)
        temp = zeros(N, N);
        temp(:, :) = pre_compute_values(ii, :, :);
        res(ii) = norm(R - temp, 2);
    else
        res(ii) = Inf;
    end
end
end
function [K] = compute_K_sur(X, L)
[N, D] = size(X);
K = zeros(N, N);
for ii = 1 : N;
    nl = length(L(L == L(ii)));
    for jj = 1 : N
        if L(ii) == L(jj)
            K(ii, jj) = 1/nl;
        end
    end
end
end
function [K] = compute_K_ulsur(X)
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(7+1, :);
sita = sita' *  sita;
K = exp(-(distance.^2)./sita);
end