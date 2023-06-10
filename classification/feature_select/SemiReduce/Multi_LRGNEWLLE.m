function [ Ws, V ] = Multi_LRGNEWLLE( XL, L, XUL, num, para)
%LRGNEWLLE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    para.a = 0.2;
    para.b = 0.2;
    num = 3;
end
if nargin == 4
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
%X =[XL;XUL];

%参数设置
K = 5;
% a = 1.1;
% b = 0.8;
% a = para.a;
% b = para.b;
a = 0.1;
b = 0.1;
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
LV = XL * LB* XL'+ b.* X * Lf * X';
RV = XL * LW* XL' + a.* X * M * X';
Ws = cell(num +1, 1);
WiWiT = 0;
aa =  trace(XL * LW* XL');
W_A = [];
V_A = [];
for ii = 1 :num
a1 = 0;
if WiWiT ~= 0
   a1 = aa /trace(WiWiT) * 0.3; 
end
RV1 =  RV +   WiWiT  * a1;
% LV = XL * LB* XL'+ b.* X * Lf * X';
% RV =  XL * LW* XL' + a.* X * M * X';
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig(LV, RV1);
% [~,D,V] = svd(pinv(RV1)*LV);
[sorted,index] = sort(diag(D),'descend');
W=V;
W_A = [W_A, W];
V_A = [V_A, diag(D)'];
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
WiWiT = WiWiT + W*W';
%WiWiT = WiWiT + W(:, 10:50)*W(:, 10:50)';

V = diag(D);
% if ii > 1
%     W1 =  Ws{1};
%     W(:,11:end) = W(:, 1:end-10);
%     W(:, 1:10) = W1(:, 1 : 10);
% end
 Ws{ii} = W;
 
end
[sorted,index] = sort(V_A,'descend');
W_A = W_A(:, index);
V_A = V_A(index);
Ws{num + 1} = W_A;
end

% function [ GF, GN ] = get_graphnf( X, K )
% %GET_GRAPHNF Summary of this function goes here
% %   Detailed explanation goes here
% [D,N] = size(X);
% %fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
% % STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
% %fprintf(1,'-->Finding %d nearest neighbours.\n',K);
% 
% X2 = sum(X.^2,1);
% distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
% 
% [sorted,index] = sort(distance);
% GN=eye(N,N);%邻接矩阵
% %S=sparse(S);
% GF=(ones(N,N)-eye(N,N))/(N-K);%非邻接矩阵
% %A=sparse(A);
% for i=1:N
%     for j=2:1+K
%         GN(i,index(i,j))=1;
%         GN(index(i,j),i)=1;
%         GF(i,index(i,j))=0;
%         GF(index(i,j),i)=0;
%     end
%     GN(i, :) = GN(i, :)./sum(GN(i, :));
%     GF(i, :) = GF(i, :)./sum(GF(i, :));
% end
% GN=GN-diag(diag(GN));
% GF=GF-diag(diag(GF));
% end
