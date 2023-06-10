function [W ] = LSPD( XL, L, XUL )
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

%参数设置
K = 5;
a = 0.1;
b = 0.1;
X = X';
XL = XL';
%第一步：计算M
[M] = getLLEGRAPH(X,K);
%M = M.^2;
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_graphnf(X, 5);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;

[SW, SB] = LSPF(XL', L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
% LV = b.* X * Lf * X' +  XL * LB * XL';
% RV =  a.* X * M * X' +  XL * LW * XL'; 

LB =  XL * LB * XL';
LW = XL * LW * XL';
Lf = X * Lf * X' ;
M = X * M * X';

LV = LB/trace(LB) + b.* Lf/trace(Lf);
RV =  LW/trace(LW) + b * M /trace(M );

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

function [SW, SB] = LSPF(XL, L)
K = 3;
%K = 2 * K;
distance = pdist2(XL, XL);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
t = mean(mean(sita));
P = exp(-distance.^2./(t));
[N, D] = size(XL);
G = zeros(N, N);
for ii = 1 : N
    for jj = 1 : N
        if L(ii) == L(jj)
            G(ii, jj) = P(ii, jj) * exp(P(ii, jj) + 1);
        else
            G(ii, jj) = -P(ii, jj) * exp(P(ii, jj) - 1);
        end
    end
end

SW = zeros(N, N);
SB = zeros(N, N);
for ii = 1 : N
    NSW = G(ii, :);
    NSW = NSW(L == L(ii));
    inx = 1:length(L);
    inx = inx(L == L(ii));
    [sorted,index] = sort(NSW, 'descend');
    neighborhood = inx(index(1:K));
    SW(neighborhood, ii) = SW(neighborhood, ii) + G(neighborhood, ii);
    SW(ii, neighborhood) = SW(ii, neighborhood) + G(ii, neighborhood);
    %SW(ii, :) = SW(ii, :)./sum(SW(ii, :));
    
    NSB = G(ii, :);
    NSB = NSB(L ~= L(ii));
    inx = 1:length(L);
    inx = inx(L ~= L(ii));
    [sorted,index] = sort(NSB, 'descend');
    neighborhood = inx(index(1:K * 3));
    SB(neighborhood, ii) = SB(neighborhood, ii) - G(neighborhood, ii);
    SB(ii, neighborhood) = SB(ii, neighborhood) - G(ii, neighborhood);   
    %SB(ii, :) = SB(ii, :)./sum(SB(ii, :));
end
% SW = SW + SW'- SW' * SW;
% SB = SB + SB'- SB' * SB;
% SW=SW-diag(diag(SW));
% SB=SB-diag(diag(SB));
end