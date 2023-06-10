function [ W ] = OptimizationL21Common( X, Y )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
X = X';
gama = 0.1;
n = size(X, 2);
I = eye(n, n);
A = [X' gama * I ];
m = size(A, 2);
Dt = eye(m, m);
for t=1:15
    C = (A * Dt * A')^-1;
    Ut = Dt * A' * C * Y;
    Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
 W = Ut(1:size(X, 1), :);
end

