% Code for SSDR
% always alpha=1,beta=20
function [newData,transMatrix] = SSDR(data,dimr,alpha,beta)
% Input
%      - data
%              data.X  : Origin data, every row denotes a data vector
%              data.ML : must-link constraints dataset
%              data.CL : cannot-link constraints dataset
%      - dimr  : The embedded dimension
%      - alpha : The parameter, default=1
%      - beta  : The parameter, default=20
% Output
%      - newData     : The embedded data
%      - transMatrix : The transfer matrix W
%
% Reference: D. Zhang, Z.-H. Zhou, and S. Chen. Semi-supervised dimensionality reduction. 
%            In: Proceedings of the 7th SIAM International Conference on Data Mining (SDM'07), 
%            Minneapolis, MN, 2007, pp.629-634.
%
% ATTN:  This package is free for academic usage. You can run it at your own risk.
%        For other purposes, please contact Dr. Daoqiang Zhang (dqzhang@nuaa.edu.cn).

[n,d] = size(data.X);

if ~exist('dimr','var') | dimr<1
    dimr = 1;
end
if ~exist('alpha','var')
    alpha = 1;
end
if ~exist('beta','var')
    beta = 20;
end

X = data.X;
ML = data.ML;
CL = data.CL;

nM = size(ML,1);
nC = size(CL,1);

% ======= Computing Graph Laplacian Matrix
S = (1/(n*n))*ones(n,n);
for i=1:nM
    tx = ML(i,1);
    ty = ML(i,2);
    S(tx,ty) = S(tx,ty) - beta/nM;
    S(ty,tx) = S(ty,tx) - beta/nM;
end 
for i=1:nC
    tx = CL(i,1);
    ty = CL(i,2);
    S(tx,ty) = S(tx,ty) + alpha/nC;
    S(ty,tx) = S(ty,tx) + alpha/nC;
end 
clear i tx ty nM nC ML CL;

S = max(S,S');% make sure S is symmetrical
D = sum(S,2);
D = diag(D);
L =  D-S;
L(isnan(L)) = 0; D(isnan(D)) = 0; 
L(isinf(L)) = 0; D(isinf(D)) = 0;
if d > 500 % big dimension
    [X,junk] = PCA(X,dimr,0.98);% Preprocess with PCA, PCAradio=0.98
    d= size(X,2);
end
LPrime = X'*L*X;

% ========= Computing eigen problem, max: W'*(X'*L*X)*W, st. W'*W=I
I = eye(d,d);
[eigvec, eigval] = eig( LPrime,I);
eigval = diag(eigval);
[eigval, index] = sort(eigval, 'descend');
eigval = eigval(1:dimr);
transMatrix = eigvec(:,index(1:dimr));% n x dimr matrix
newData = X * transMatrix;

end

