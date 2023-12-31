function [ Ws, V ] = Multi_LRGPSSDR_GCMCTRACE(  XL, L, XUL, num  )
%LOCATELRGPSSDR Summary of this function goes here
%   Detailed explanation goes here
%特征归一化
if nargin < 4
    num = 3;
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

%参数设置
UL = unique(L);
K = 10000;
for z = 1 : length(UL)
    if ( K > length(L(L == UL(z))))
         K = length(L(L == UL(z)));
    end
end
a = 0.3;
b = 0.3;
X = X';
XL = XL';


%第一步：计算M
%[M] = getLLEGRAPH(X,K);
[ M ] = get_GCMC_W(  X', K ); 
%M = M.^2;
%第二步：计算k近远拉普拉斯图
 [GF, GN] = get_graphnf(X, K);
%[ GF ] = get_GCMC_W(  X', 10 ); 
 Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%第三步：计算类内类间散布矩阵
[SB, SW, ST] = get_weight_SW_SB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;

LB = XL * LB* XL';
Lf = X * Lf * X';
LW = XL * LW* XL';
M = X * M * X';
 % M1 = X * M1 * X';
LV = LB/trace(LB) + b.* Lf/trace(Lf)';
RV = LW/trace(LW) + a.* M/trace(M);
Ws = cell(num, 1);
 WiWiT = 0;
for ii = 1 :num 

RV1  = RV +  ii * WiWiT /(trace(WiWiT ) + 0.0000001);% + a.* M1/trace(M1);
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig((LV), (RV1));
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
WiWiT = WiWiT + W*W';
% for i = 1:size(W,2)
%     W(:,i) = W(:,i)./norm(W(:,i));
% end
V = diag(D);
 Ws{ii} = W;
end
end


