function [ W ] = OptimizationL21L22Common_co1( X, Y, Y1, Para)
%OPTIMIZATIONL21L22COMMON Summary of this function goes here
%   Detailed explanation goes here
miu = 0.1;
lamda = Para.lambda;
beta = Para.beta;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:特征选择时， 样本个数*样本维数
%X稀疏表示时，样本维数*样本个数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut1 = eye(m, m);
Ut2 = eye(m, m);
Ut3 = eye(m, m);
XtY = X'* Y;
XtY1 = X'* Y1;
XtX = X'* X ;
for t=1:15
    
    
    W1 = (Ut1*Ut3 * XtX + lamda * Ut1  + beta * Ut3+ 0.00000001 * I1)^-1 * Ut1*Ut3 * XtY1;
    Ut1 = diag(sqrt(sum(W1.^2, 2)).*2);
%     if t== 1
%         Ut3 = Ut1;
%     else
%         Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
%     end
    W2 = (Ut2*Ut3 * XtX+ lamda * Ut2 + beta * Ut3+ 0.00000001 * I1)^-1 * Ut2*Ut3 * XtY;
    Ut2 = diag(sqrt(sum(W2.^2, 2)).*2);    
    W = [W1, W2];    
    Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
%     W = (Ut* XtX+ lamda * I1)^-1 * Ut * XtY;
%     Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
end

