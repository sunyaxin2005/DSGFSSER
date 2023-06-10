function [ W ] = OptimizationL21EbU_Right_manifold( X, Y, L, Para )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
mean_x = mean(mean( eye(size(X, 1), size(X, 1))));
mean_L = mean(mean(abs(L)));
miu = Para.beta * mean_x/mean_L;
lamda = Para.lambda;
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
XtY = X'* Y;
XtCX = X'* (I + miu * L) * X;
Ut = eye(m, m);
AB  = (Ut* XtCX+ lamda * I1+ 0.00000001 * I1)^-1 * Ut * XtY;
LV = X'*Y*Y'*X;
I2 = eye(size(Y, 2), size(Y, 2));
for t=1:15
    RV = XtCX + lamda * Ut;
    [V,D]=eig(LV, RV);
    [sorted,index] = sort(diag(D),'descend');
    A = zeros(m, size(Y, 2));
    for i=1:size(Y, 2)
        A(:,i) = V(:,index(i));
    end
    %B = 
    B = (A'*Ut*XtCX*A+lamda*(A'*A)+ 0.00000001 * I2)^-1 * (A'*Ut*XtY);
    AB = A*B;
    Ut = diag(sqrt(sum(AB.^2, 2)).*2);  
end
W = AB;
end