function [ W  ] = OptimizationL21Manifold(  X, Y, L, Para)
%OPTIMIZATIONL21MANIFOLD Summary of this function goes here
%   Detailed explanation goes here
miu = Para.beta;
lamda = Para.lambda;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut = eye(m, m);
XtY =  X'* Y;
XtCX = X'* (I + miu * L) * X;
for t=1:15
    W = (Ut* XtCX + lamda * I1)^-1 * Ut * XtY;
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end

end

