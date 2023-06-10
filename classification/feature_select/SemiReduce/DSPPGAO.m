function [ W ] = DSPPGAO(  XL, L, XUL  )
%DSPPGAO Summary oaf this function goes here
%   Detailed explanation goes here
%Discriminative sparsity preserving projections for image recognition (Pattern Recognition 2015)
X =[XL;XUL];
[ S ] = get_DSPNE(  XL, L );
X1 = XL';
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
% %第一步得到稀疏S
% [ S ] = get_DSPNE(  XL, L );
%第三步：计算类内类间散布矩阵
% [ sb, sw ] = get_sbsw( XL, L );
XL = XL';
% 
[GB] = get_B(XL, L, 6);
DB = sum(GB(:,:),2);
LB = diag(DB)-GB;
[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
LV = XL * LB* XL';
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


function [GN] = get_B(X, L, K)
%b = 1 if x_j belong N(K) and L(xi) ~= L(xj); b=0, otherwise
[D,N] = size(X);
%fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
%fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
GN=zeros(N,N);%邻接矩阵
%S=sparse(S);
%A=sparse(A);
for i=1:N
    for j=2:1+K
        if L(i) ~= L(index(i,j))
            GN(i,index(i,j))=1;
            GN(index(i,j),i)=1;
        end
    end
end
end

