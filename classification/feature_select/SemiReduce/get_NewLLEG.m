function [ L ] = get_NewLLEG( X, K )
%NEWLLEG Summary of this function goes here
%   Detailed explanation goes here
%定义U，按离第rank_index个元素的远近进行排序，那么第rank_index对角的值赋值非常大，其它为1
N = size(X, 1);
X = X';
lambuda = 1000;
a = 0.1;
gamma = 0.5;
distance = pdist2(X', X');
[sorted,index] = sort(distance);
neighborhood = index(:, 1:(1+K));
H = zeros(K + 1, 1);
H = H + 1;
I = eye(K+1);
H = I - H * H';

S = [];
for ii = 1 : N
    Si = zeros(N, K + 1); 
    for jj = 1 : K + 1
        Si(neighborhood(ii, jj), jj) = 1;
    end
    S = [S, Si];
end

A = zeros(N * (K + 1),  N * (K + 1));
for ii = 1 : N
    Li = lambuda * H *(H * X(:, neighborhood(ii, :))' * X(:, neighborhood(ii, :)) * H + lambuda * I)^(-1) * H;
    start_i = (ii - 1) * (K + 1) + 1;
    end_i = start_i + (K + 1) - 1;
    A(start_i:end_i,  start_i:end_i) = Li;
end

L = S * A * S';

for ii = 1 : N
    L(ii, :) = L(ii, :)./sum( L(ii, :));
end
L = L + L' - L * L';
end

