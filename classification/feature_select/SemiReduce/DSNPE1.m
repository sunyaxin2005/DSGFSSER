function [ W ] = DSNPE1(  XL, L, XUL)
%DSNPE Summary of this function goes here
%   Detailed explanation goes here
X =[XL;XUL];
[ S ] = get_DSPNE1(  X );
X1 = X';
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;
%参数设置
X = X';
beta = 0.1;
a = 0.1;
b = 0.1;
%第一步得到稀疏S
% [ S ] = get_DSPNE1(  X );
%第三步：计算类内类间散布矩阵
% [ sb, sw ] = get_sbsw( XL, L );
XL = XL';
% 
[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
[GF, GN] = get_graphnf(X, 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
LV = XL * LB* XL' + b.* X * Lf * X';
RV =  XL * LW* XL' + a.* X1 * S' *X1';
%  LV = sb;
%  RV = sw;
%LV = XL * S * XL'- beta * (sb- sw);
% RV = XL * XL';
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig(LV, RV);
diagD = diag(D);
diagD(isinf(diagD)) = 0;
[sorted,index] = sort(diagD,'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end

for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
end

