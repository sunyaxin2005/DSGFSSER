function [ W ] = OptimizationL21EbU_Right( X, Y, Para )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:

lamda = Para.lambda;
n = size(X, 1);
m = size(X, 2);
I1 = eye(m, m);
XtY = X'* Y;
XtX = X'* X ;
Ut = eye(m, m);
AB  = (Ut* XtX+ lamda * I1 + 0.00000001 * I1)^-1 * Ut * XtY;
LV = X'*Y*Y'*X;
I2 = eye(size(Y, 2), size(Y, 2));
for t=1:15
    RV = XtX + lamda * Ut;
    [V,D]=eig(LV, RV);
    [sorted,index] = sort(diag(D),'descend');
    A = zeros(m, size(Y, 2));
    for i=1:size(Y, 2)
        A(:,i) = V(:,index(i));
    end
    %B = 
    B = (A'*Ut*XtX*A+lamda*(A'*A) + 0.00000001 * I2)^-1 * (A'*Ut*XtY);
    AB = A*B;
    Ut = diag(sqrt(sum(AB.^2, 2)).*2);  
end
W = AB;
end