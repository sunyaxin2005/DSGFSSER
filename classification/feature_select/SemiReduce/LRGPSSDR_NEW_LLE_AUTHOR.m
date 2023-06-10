function [ W  ] = LRGPSSDR_NEW_LLE_AUTHOR(  XL, L, XUL, para )
%LRGPSSDR_NEW_LLE_AUTHOR Summary of this function goes here
%   Detailed explanation goes here
if nargin <= 3
   k = 5;
   a = 0.1;
else
   k = para.k;
   a = para.a;
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

% [ LA ] = get_NewLLEG( X, k);
% %第二步：计算k近远拉普拉斯图
% [GF, GN] = get_graphnf(X, 5);
% Df = sum(GF(:,:),2);
% Lf = diag(Df)-GF;
% Dn = sum(GN(:,:),2);
% Ln = diag(Dn)-GN;
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
%[SB, SW, ST] = get_SWSB(XL, L);
% [ SB, SW, ST ] = get_new_SW_SB_UNNORMAL( XL, L );
% DB = sum(SB(:,:),2);
% LB = diag(DB)-SB;
% DW = sum(SW(:,:),2);
% LW = diag(DW)-SW;
% DT = sum(ST(:,:),2);
% LT = diag(DT)-ST;
LB = XL * LB* XL';
LW = XL * LW* XL';
LA = X*M*X';
% Lf = X*Lf*X';

RV = LW/trace(LW) + a * LA/trace(LA) ;
SB = LB/trace(LB);% + b * Lf/trace(Lf) ;
LV = SB;

[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end

for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
end

