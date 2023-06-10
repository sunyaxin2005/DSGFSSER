function [ W ] = OptimizationL21L22Common_co3( X, Y, Para)
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
miu = 0.1;
lamda = Para.lambda;
beta = Para.beta;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:特征选择时， 样本个数*样本维数
%X稀疏表示时，样本维数*样本个数
%Y:
%X = X';
Y1 = X;
% for i = 1:size(Y1,2)
%     Y1(:,i) = Y1(:,i)./norm(Y1(:,i));
% end
% for i = 1:size(Y,2)
%     Y(:,i) = Y(:,i)./norm(Y(:,i));
% end
% for i = 1:size(X,2)
%     X(:,i) = X(:,i)./norm(X(:,i));
% end
b = (size(Y, 2)/size(X, 2))/10;
%Y = repmat(Y, 1,b);
%Y = Y(:, 1: size(X, 2));
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut1 = eye(m, m);
Ut2 = eye(m, m);
Ut3 = eye(m, m);
XtY = X'* Y;
XtY1 = X'* Y1;
a = diag(diag(XtY1));
XtY1 = XtY1 - a;
XtX = X'* X ;
for t=1:15
    W1 = (b*Ut1*Ut3 * XtX+ b*lamda * Ut1 + beta * Ut3 + 0.00000001 * I1)^-1 * Ut1*Ut3 * XtY1;
    Ut1 = diag(sqrt(sum(W1.^2, 2)).*2);
%     if t== 1
%         Ut3 = Ut1;
%     else
%         Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
%     end
    W2 = (Ut2*Ut3 * XtX+ lamda * Ut2 + beta * Ut3+ 0.00000001 * I1)^-1 * Ut2*Ut3 * XtY;
    Ut2 = diag(sqrt(sum(W2.^2, 2)).*2);    
    W = [b*W1, W2];    
    Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
%     W = (Ut* XtX+ lamda * I1)^-1 * Ut * XtY;
%     Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
W = W2;
end

