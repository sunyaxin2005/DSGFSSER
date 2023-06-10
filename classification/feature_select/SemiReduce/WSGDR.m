function [ W ] = WSGDR(  XL, L, XUL )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if(size(XL, 2) == 1024)
    [XL] = vector2symvector( XL );
    [XL] = vector2symvector( XL );
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

%参数设置
K = 1;
a = 0.1;
b = 0.1;
X = X';
XL = XL';
%第一步得到稀疏S
%[ S1 ] = get_DSPNE(  XL', L );
[ S2 ] = get_NewLLEG( X',K );
%[ S ] = get_GCMC_W(  X', K ); 
%第三步：计算类内类间散布矩阵
% [ sb, sw ] = get_sbsw( XL, L );
%
LG = zeros(size(X, 1), 1) + 1;
LG(1:floor(size(X, 1)/2)) = 1;
LG(size(X, 1):-1:size(X, 1)-floor(size(X, 1)/2) + 1) = -1;
LG = LG*LG';

[GF, GN] = get_graphnf(X, 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;

[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
LB = LB/(trace(LB));

LB = XL * LB* XL';
LW =  XL * LW* XL';
Lf = X * Lf * X';
%LS =  X * S'*X';
LS2 =  X * S2'*X';
%LS1 =  XL * S1'*XL';

LW = LW/(trace(LW));
LB = LB/(trace(LB));
Lf = Lf/trace(Lf);
%LS = LS/trace(LS);
%LS1 = LS1/trace(LS1);
LS2 = LS2/trace(LS2);
LG = LG/trace(LG);
LV = LB + b.* Lf;
RV =  LW + a.* LG + a.*LS2;
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

