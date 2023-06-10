function [ GF, GN ] = get_graphnf( X, K )
%GET_GRAPHNF Summary of this function goes here
%   Detailed explanation goes here
[D,N] = size(X);
%fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
%fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
GN=eye(N,N);%ÁÚ½Ó¾ØÕó
%S=sparse(S);
GF=(ones(N,N)-eye(N,N))/(N-K);%·ÇÁÚ½Ó¾ØÕó
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

