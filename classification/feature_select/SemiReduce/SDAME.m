function [ W ] = SDAME(  XL, L, XUL  )
%LOCATELRGPSSDR Summary of this function goes here
%   Detailed explanation goes here
%特征归一化
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
UL = unique(L);
K = 10000;
for z = 1 : length(UL)
    if ( K > length(L(L == UL(z))))
         K = length(L(L == UL(z)));
    end
end
a = 0.1;
b = 0.1;
X = X';
XL = XL';
%第一步：计算M
%[M] = getLLEGRAPH(X,K);
[ M ] = get_LE_W(  X', K );
%M = M.^2;
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_graphnf(X, 5);
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
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
LV = XL * LB* XL';% + b.* X * Lf * X';
RV =  XL * LW* XL' + a.* X * M * X';
%[U1,D1,QP] = svd(inv(RV)*LV);
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
function [ M ] = get_LE_W(  X, K )
%GET_GCMC_W Summary of this function goes here
%   Detailed explanation goes here
distance = pdist2(X, X);
sita = mean(mean(distance));
W = exp(-distance.^2/(sita * sita));
%W = distance .* W;
for ii=1:size(W, 2)
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;
M = (W+W')/2;

end

