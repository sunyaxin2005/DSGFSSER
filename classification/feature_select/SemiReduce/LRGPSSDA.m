function [ W ] = LRGPSSDA( XL, L, XUL )
%LRGPSSDR Summary of this function goes here
%   Detailed explanation goes here
%X:points * dimension
%L:数据的Lable
K = 5;
a = 0.1;
b= 0.1;
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
%第一步：计算M
[M] = getLLEGRAPH(X,K);
%第二步：计算k近远拉普拉斯图
[GF, GN] = get_gf(X, K);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%第三步：计算类内类间散布矩阵
[Sb, Sw] = compute_swsb(XL, L);
LV = Sb + b.* X * Lf * X';
RV = Sw + a.* X * M * X';
%[U1,D1,QP] = svd(inv(RV)*LV);
%[V,D]=eig(RV, LV);
[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
%W = normc(U1(:,1:end));
%W = D;
end

%计算近远邻图
function [GF, GN] = get_gf(X, K)
[D,N] = size(X);
%fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
%fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
GN=eye(N,N);%邻接矩阵
%S=sparse(S);
GF=ones(N,N)-eye(N,N);%非邻接矩阵
%A=sparse(A);
for i=1:N
    for j=2:1+K
        GN(i,index(i,j))=1;
        GN(index(i,j),i)=1;
        GF(i,index(i,j))=0;
        GF(index(i,j),i)=0;
    end
end
end

function [Sb, Sw] = compute_swsb(X, labels)
X = X';
%fprintf(1,'compute_swsb running on %d points in %d dimensions\n',N,D);
[classes, bar, labels] = unique(labels);
nc = length(classes);
Sw = zeros(size(X, 2), size(X, 2));
    
% Compute total covariance matrix
St = cov(X);

	% Sum over classes
for i=1:nc
        
    % Get all instances with class i
    cur_X = X(labels == i,:);

	% Update within-class scatter
	C = cov(cur_X);
	p = size(cur_X, 1) / (length(labels) - 1);
	Sw = Sw + (p * C);
end
    
    % Compute between class scatter
Sb       = St - Sw;
Sb(isnan(Sb)) = 0; Sw(isnan(Sw)) = 0;
Sb(isinf(Sb)) = 0; Sw(isinf(Sw)) = 0;
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