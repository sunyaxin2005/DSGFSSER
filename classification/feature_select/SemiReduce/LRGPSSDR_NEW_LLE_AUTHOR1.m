function [ W  ] = LRGPSSDR_NEW_LLE_AUTHOR1(  XL, L, XUL , para)
%LRGPSSDR_NEW_LLE_AUTHOR Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    para.a = 0.2;
    para.b = 0.2;
else
 %   k = para.k;
    a = para.a;
    b = para.b;
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

X = X';
XL = XL';

%第一步：计算M
[M] = get_NewLLEG(X',2);
%M = M.^2;
%第二步：计算k近远拉普拉斯图
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_graphnf(X, 5);
%[ SF ] = get_graph_farL(  XL, L,k);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
%第三步：计算类内类间散布矩阵
[SB, SW, ST] = get_SWSB(XL, L);
%[SB, SW, ST] = get_new_SW_SB_UNNORMAL(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;

LV = XL * LB* XL'/trace(XL * LB* XL');%+ b.* X * Lf * X'/trace(X * Lf * X');
RV =  XL * LW* XL'/trace(XL * LW* XL') + a.* X * M * X'/trace(X * M * X');
% LV = XL * LB* XL';%+ b.* X * Lf * X';
% RV =  XL * LW* XL' + a.* X * M * X';
%[U1,D1,QP] = svd(inv(RV)*LV);
%[V,D] = eig(LV, RV);
[V,D]=mySVD(inv(RV)*LV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end

for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end

end

