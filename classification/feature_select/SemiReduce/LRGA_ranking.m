% Ranking with Local Regression and Global Alignment
function [F] = LRGA_ranking(para, X, Y)
% para: parameters of the LRGA ranking algorithm
%   para.k:     neighborhood number to construct the graph, this parameter
%               should be a small number, such as 5, 10 and 15.
%   para.lamda: regularization parameter in Eq.(10) to control the capacity
%               of the local regression models. In most cases, the ranking proformance is not sensitive to this parameter. However, this parameter
%               should be a large constant (e.g. 10000) to make the performance stable.
% X: each column is a data point
% Y: a vector of size n where n is the number of data. If X_i is a query, Y_i=1; else Y_i=0;
% F: the output ranking values.

% Reference: 
% Yi Yang, Dong Xu, Feiping Nie, Jiebo Luo and Yueting Zhuang. Ranking with Local Regression and Global Alignment for Cross-Media Retrieval. ACM MM 2009



L = Laplacian_LRGA(X, para); % Computing the Laplacian matrix. The Laplacian matrix can be computed off line to speed up the retrieval.

n = size(X,2);
label_idx = find(Y == 1);
unlabel_U = 1;
label_U = 10000;
U = unlabel_U*ones(n,1); U(label_idx) = label_U; U = diag(U);
M = U + L;
F = M\(U*Y);



function L = Laplacian_LRGA(X, para)
k = para.k;
lamda = para.lamda;

[Dim, n] = size(X);
Lc = eye(k) - 1/k*ones(k);
A = spalloc(n*k,n*k,5*n*k);
S = spalloc(n,n*k,5*n*k);
for i = 1:n
    dis = repmat(X(:,i),1,n) - X;
    dis = sum(dis.*dis);
    [dumb, nnidx] = sort(dis);
    Xi = X(:,nnidx(1:k));
    Xi = Xi*Lc;
    if Dim > k
        Ai = inv(lamda*eye(k) + Xi'*Xi);
        Ai = lamda*Lc*Ai*Lc;
    else
        Ai = Lc - Xi'*inv(lamda*eye(Dim) + Xi*Xi')*Xi;
    end;
    lidx = (i-1)*k+1:(i-1)*k+k;
    A(lidx, lidx) = Ai;
    S(nnidx(1:k),lidx) = eye(k);
end;
L = lamda*S*A*S';
    
    

