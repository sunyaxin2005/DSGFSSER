function [ W ] = OptimizationL21Eb_Right( X, Y, Para )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:��������*����ά��
%Y:

lamda = Para.lambda;
n = size(X, 1);
m = size(X, 2);
e = zeros(n, 1)+1;
I1 = eye(m, m);
etY = e'* Y;
etX = e'* X;
XtY = X'* Y;
XtX = X'* X ;
Ut = eye(m, m);
AB  = (Ut* XtX+ lamda * I1+ 0.00000001 * I1)^-1 * Ut * XtY;
I2 = eye(size(Y, 2), size(Y, 2));
for t=1:15
    b = (etY-etX*AB)./n;
    LV = X'*(Y-e*b)*(Y-e*b)'*X;
    RV = XtX + lamda * Ut;
    [V,D]=eig(LV, RV);
    [sorted,index] = sort(diag(D),'descend');
    A = zeros(m, size(Y, 2));
    for i=1:size(Y, 2)
        A(:,i) = V(:,index(i));
    end
    %B = 
    B = (A'*Ut*XtX*A+lamda*(A'*A)+ 0.00000001 * I2)^-1 * (A'*Ut*XtY-(A'*Ut*X'*e*b));
    AB = A*B;
    Ut = diag(sqrt(sum(AB.^2, 2)).*2);  
end
W = AB;
end

