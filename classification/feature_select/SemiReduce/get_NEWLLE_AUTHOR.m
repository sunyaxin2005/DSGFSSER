function [ L ] = get_NEWLLE_AUTHOR( X, K  )
%GET_NEWLLE_AUTHOR Summary of this function goes here
%   Detailed explanation goes here
%数据每一列为一个样本
k = K;
lamda = 1000;

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
L = S*A*S';
%L = L + L' - L * L';
end

% k = K;
% lamda = 1000;
% 
% [Dim, n] = size(X);
% Lc = eye(k) - 1/k*ones(k);
% A = spalloc(n*k,n*k,5*n*k);
% S = spalloc(n,n*k,5*n*k);
% for i = 1:n
%     dis = repmat(X(:,i),1,n) - X;
%     dis = sum(dis.*dis);
%     [dumb, nnidx] = sort(dis);
%     Xi = X(:,nnidx(1:k));
%     Xi = Xi*Lc;
%     if Dim > k
%         Ai = inv(lamda*eye(k) + Xi'*Xi);
%         Ai = lamda*Lc*Ai*Lc;
%     else
%         Ai = Lc - Xi'*inv(lamda*eye(Dim) + Xi*Xi')*Xi;
%     end;
%     lidx = (i-1)*k+1:(i-1)*k+k;
%     A(lidx, lidx) = Ai;
%     S(nnidx(1:k),lidx) = eye(k);
% end;
% L = S*A*S';

