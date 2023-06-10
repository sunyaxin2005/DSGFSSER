function [ W ] = OptimizationL21Eb( X, Y, Para )
%OPTIMIZATIONL21COMMON Summary of this function goes here
%   Detailed explanation goes here
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:

lamda = Para.lambda;

n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
e = zeros(n, 1)+1;
I1 = eye(m, m);
etY = e'* Y;
etX = e'* X;
XtY = X'* Y;
XtX = X'* X ;
Ut = eye(m, m);
AB  = (Ut* XtX+ lamda * I1)^-1 * Ut * XtY;
H = I -e*e'/n;
H = H*H';
XtHHYYtHHX = X' *H *H*(Y*Y')*H*H*X;
XtHHX = X'*H*H*X;
for t=1:15
    b = (etY-etX*AB)./n;
    LV = XtHHYYtHHX;
    RV = XtHHX + lamda * Ut;
    [V,D]=eig(LV, RV);
    [sorted,index] = sort(diag(D),'descend');
    A = zeros(m, size(Y, 2));
    for i=1:size(Y, 2)
        A(:,i) = V(:,index(i));
    end
    %B = 
    B = (A'*Ut*XtX*A+lamda*(A'*A))^-1 * (A'*Ut*XtY-(A'*Ut*X'*e*b));
    AB = A*B;
    Ut = diag(sqrt(sum(AB.^2, 2)).*2);  
end
W = AB;
end

