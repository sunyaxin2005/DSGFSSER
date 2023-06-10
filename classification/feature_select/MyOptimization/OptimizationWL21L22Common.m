function [ W ] = OptimizationWL21L22Common( X, Y, S, Para)
%OPTIMIZATIONL21L22COMMON Summary of this function goes here
%   Detailed explanation goes here
miu = 0.1;
lamda = Para.lambda;
%min_W,E sum_i=1,sizeY, ||X'W-Yi||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut = eye(m, m);
XtY = X'* Y;
XtX = X'* X ;
for ii = 1 : size(S, 2)
   XtY(:, ii) = XtY(:, ii) * S(ii); 
end
XtY = XtY * size(S, 2);
for t=1:15
    W = (Ut* XtX * sum(S)+ lamda * I1)^-1 * Ut * XtY;
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
end

