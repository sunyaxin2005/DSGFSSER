function [ W ] = SSLFDA( XL, L,  XUL)
%SSLFDA Summary of this function goes here
%   Detailed explanation goes here
%ESELF Summary of this function goes here
%   Detailed explanation goes here
X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./stdd;
K = 1;
a = 0.1;
b = 0.1;
X = X';
XL = XL';
XUL = XUL';
%第一步：计算M
[M] = least_squres(X,K);
%M = M.^2;
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_gf(X, 10);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%第三步：计算类内类间散布矩阵
[SB, SW, ST] = compute_swsb(XL, L, XUL);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
X = [XL, XUL];
LV = X * LB* X'- X * LW * X';
St = cov(X');
RV =  St;
%[U1,D1,QP] = svd(inv(RV)*LV);

[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
%W = U1;
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
end

%计算近远邻图
%计算近远邻图
function [GF, GN] = get_gf(X, K)
[D,N] = size(X);
%fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
%fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
GN=eye(N,N);%邻接矩阵
%S=sparse(S);
GF=(ones(N,N)-eye(N,N))/(N-K);%非邻接矩阵
%A=sparse(A);
for i=1:N
    for j=2:1+K
        GN(i,index(i,j))=1/K;
        GN(index(i,j),i)=1/K;
        GF(i,index(i,j))=0;
        GF(index(i,j),i)=0;
    end
end
end

function [SB, SW, ST] = compute_swsb(XL, L, XUL)
X = [XL,XUL];
K = 5;
X = X';
[N, D] = size(X);
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
A = exp(-(distance.^2)./sita);
diagA = diag(ones(size(distance, 2), 1));
A = A - diagA;
neighborhood = index(2:(1+K),:);
% A1 = A;
% for ii = 1 : size(neighborhood, 1)
%     A1(ii, neighborhood(ii, :)) = 0; 
%     A1(neighborhood(ii, :), ii) = 0; 
% end
% A1 = A - A1;
SW = zeros(N, N);
UL = unique(L);
nClass = length(UL);
%赋值有label的
for ii = 1 : nClass
    SW(L == UL(ii), L == UL(ii)) = A(L == UL(ii), L == UL(ii))/length(L(L == UL(ii)));
end% %赋值没有label的
%赋值没有label的
for ii = 1 : N
    SW(neighborhood(:, ii), ii) = SW(neighborhood(:, ii), ii) + A(neighborhood(:, ii), ii)/N;
    SW(ii, neighborhood(:, ii)) = SW(ii, neighborhood(:, ii)) + A(ii, neighborhood(:, ii))/N;
end
SW=SW-diag(diag(SW));
%赋值有label的
SB = zeros(N, N);
for ii = 1 : nClass
    SB(L == UL(ii), L == UL(ii)) = -A(L == UL(ii), L == UL(ii))*(1/length(L(L == UL(ii))) - 1/N);
    SB(L == UL(ii), L ~= UL(ii)) = A(L == UL(ii), L ~= UL(ii))/N;
    SB(L ~= UL(ii), L == UL(ii)) = A(L ~= UL(ii), L == UL(ii))/N;
end% %赋值没有label的
SB=SB-diag(diag(SB));
ST = A;
end

function [M] = least_squres(X,K)
% LLE ALGORITHM (using K nearest neighbors)
%
% [Y] = lle(X,K,dmax)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% K = number of neighbors
% dmax = max embedding dimensionality
% Y = embedding as dmax x N matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,N] = size(X);
%fprintf(1,'least_squres running on %d points in %d dimensions\n',N,D);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
%fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);

% STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
%fprintf(1,'-->Solving for reconstruction weights.\n');

if(K>D) 
  %fprintf(1,'   [note: K>D; regularization will be used]\n'); 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end

W = zeros(K,N);
for ii=1:N
   z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
   C = z'*z;                                        % local covariance
   C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
   W(:,ii) = C\ones(K,1);                           % solve Cw=1
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;


% STEP 3: COMPUTE EMBEDDING FROM EIGENVECTS OF COST MATRIX M=(I-W)'(I-W)
%fprintf(1,'-->Computing embedding.\n');

% M=eye(N,N); % use a sparse matrix with storage for 4KN nonzero elements
M = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M(ii,jj) = M(ii,jj) - w';
   M(jj,ii) = M(jj,ii) - w;
   M(jj,jj) = M(jj,jj) + w*w';
end;
end
