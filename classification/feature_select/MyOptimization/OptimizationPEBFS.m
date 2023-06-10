function [ W] = OptimizationPEBFS(  X, Y, Para )
%OPTIMIZATIONPEBFS 此处显示有关此函数的摘要
%   此处显示详细说明
lamda = Para.lambda;
[n, m] = size(X);
In = eye(n, n);
l = zeros(n, 1) + 1;
X = X';
A = X*(In - l*l'./n)*X';
B = X * (In - l*l'./n) * Y;
Q = eye(m, m);
C = size(Y, 2);
for t = 1 : 15
    W = Q(:, 1 : C);
    H = Q(:, (C+1):m);
    r = trace(W' * B)/trace(W' * A * W);
    D = diag(1./(sqrt(sum(W.^2, 2)).*2 + 0.0000001));
    R = [B, r * A * H./2 + (lamda/(2 * r)) * D * H];
    [U,S,V] = svd(R);
    Q = U*V';
end


end

