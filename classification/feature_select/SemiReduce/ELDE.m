function [ W  ] = ELDE(  XL, L, XUL )
%ELDE 此处显示有关此函数的摘要
%   此处显示详细说明
%Matrix exponentialbasedsemi-superviseddiscriminantembedding for
%imageclassification Pattern recognition 2017
X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;

%参数设置
K = 1;
a = 0.1;
b = 0.1;
X = X';
XL = XL';
%第一步：计算M
[M] = getLLEGRAPH(X,K);
%M = M.^2;
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_graphnf(X, 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%第三步：计算类内类间散布矩阵
[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
% DT = sum(ST(:,:),2);
% LT = diag(DT)-ST;
% LV = XL * LB* XL' + b.* X * Lf * X';
% RV =  XL * LW* XL' + a.* X * M * X';
LB = XL * LB* XL';
Lf = X * Lf * X';
LW = XL * LW* XL';
Ln= X * Ln * X';
LV = LB/trace(LB);% + b.* Lf/trace(Lf)';
RV =  LW/trace(LW) + a.* Ln/trace(Ln);
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig(exp(LV), exp(RV));
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
V = diag(D);
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end

end

