function [ Wt ] = OptimizationL21SCM( X, Y, Para )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
lamda = Para.lambda;

X = X';
Y = Y';
[d, n] = size(X);
Ft = eye(n, n);
Dt = eye(d, d);
I1 = eye(d, d);
for t=1:15
    Wt = (Dt * X * Ft * X' + lamda * I1)^-1 * Dt * X * Ft * Y';
    Ft = sum((Wt'* X - Y).^2);
    B = sort(Ft, 'descend');
    e = B(floor(n/10));
    Ft(Ft > e) = 0;
    Ft = diag(Ft);
    Dt = diag(sqrt(sum(Wt.^2, 2)).*2);
end
end

