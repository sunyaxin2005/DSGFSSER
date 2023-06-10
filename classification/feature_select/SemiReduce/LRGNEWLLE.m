function [ W, V ] = LRGNEWLLE( XL, L, XUL, para  )
%LRGNEWLLE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    para.a = 0.2;
    para.b = 0.2;
end

 X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;

X(isnan(X)) = 0;
X(isinf(X)) = 0;

%参数设置
K = 1;
% a = 1.1;
% b = 0.8;
% a = para.a;
% b = para.b;
a = 0.3;
b = 0.3;

X = X';
XL = XL';
%第一步：计算M
[M] = get_NewLLEG(X',K);
%M = M.^2;
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_graphnf(X, 5);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%第三步：计算类内类间散布矩阵
[SB, SW, ST] = get_SWSB(XL, L);
%[SB, SW, ST] = get_new_SW_SB_UNNORMAL(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
LV = XL * LB* XL'/trace(XL * LB* XL')+ b.* X * Lf * X'/trace(X * Lf * X');
RV =  XL * LW* XL'/trace(XL * LW* XL') + a.* X * M * X'/trace(X * M * X');
% LV = XL * LB* XL'+ b.* X * Lf * X';
% RV =  XL * LW* XL' + a.* X * M * X';
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig(LV, RV);
% [~,D,V] = svd(pinv(RV)*LV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
% 
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
V = diag(D);
end

