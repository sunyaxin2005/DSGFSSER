function [ W ] = OptimizationL21L22Common( X, Y, Para)
%OPTIMIZATIONL21L22COMMON Summary of this function goes here
%   Detailed explanation goes here
miu = 0.1;
%lamda = Para.lambda;
lamda = 0.8;
%min_W,E ||X'W-Y||_F,28 + ||W||_2,1
%X:特征选择时， 样本个数*样本维数
%X稀疏表示时，样本维数*样本个数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut = eye(m, m);
XtY = X'* Y;
XtX = X'* X ;
for t=1:15
    W = (Ut* XtX+ lamda * I1)^-1 * Ut * XtY;
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
end

