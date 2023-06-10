function [ W ] = LDRCLRGPSSDE( XL, L, XUL  )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;

ata = 0.1;
XL = XL';
XUL = XUL';
C = unique(L);
C = length(C);
H = cell(C, 1);
for ii = 1 : C
    Xii = XL(:, L == ii);
    H{ii} = Xii * (Xii'*Xii)^-1 * Xii';
end

Xinter = cell(C, C);
Xintra = cell(C, 1);
for ii = 1 : C
    Xii = XL(:, L == ii);
    for jj = 1 : C
        Xinter{ii, jj} = H{jj} * Xii;
    end
    Xintra{ii} = H{ii} * Xii;
end
Eb = 0;
for ii = 1 : C
    Xii = XL(:, L == ii);
    for jj = 1 : C
        if ii == jj
            continue;
        end
        Eb = Eb+(Xii - Xinter{ii, jj}) * (Xii - Xinter{ii, jj})';
    end
end
Eb = Eb/(size(XL, 2) * (C-1));
Ew = 0;
for ii = 1 : C
    Xii = XL(:, L == ii);
    Ew = (Xii - Xintra{ii}) * (Xii - Xintra{ii})';
end
Ew = Ew/size(XL, 2) ;
I = eye(size(Ew));
LV = Eb;
RV = Ew + ata * I;
LV = LV/trace(LV);
RV = RV/trace(RV);
%参数设置
K = 1;
a = 0.5;
b = 0.5;
X = X';
%XL = XL';
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
[WSB, WSW, WST] = get_weight_SW_SB(XL, L);
WDB = sum(WSB(:,:),2);
WLB = diag(WDB)-WSB;
WDW = sum(WSW(:,:),2);
WLW = diag(WDW)-WSW;
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;

VLG = zeros(size(X, 1), 1) + 1;
LG(1:floor(size(X, 1)/2)) = 1;
LG(size(X, 1):-1:size(X, 1)-floor(size(X, 1)/2) + 1) = -1;
LG = LG*LG';
% DT = sum(ST(:,:),2);
% LT = diag(DT)-ST;
LV1 = XL * LB* XL';% + b.* X * Lf * X';
RV1 =  XL * LW* XL';% + a.* X * M * X';
LV1 = LV1/trace(LV1);
RV1 = RV1/trace(RV1);
LG = LG/trace(LG);

LV = LV + b.* LV1;
RV = RV + a.* RV1;

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

