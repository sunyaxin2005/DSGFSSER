function [W, F] = FME_unsupervise(X, L, para)
% X: each colomn is a data point
% L: Laplacian matrix, the matrix M in the paper
% para: parameters
%       para.ul: parameter associated to labeled data points
%       para.uu: parameter associated to unlabeled data points(usually be 0)
%       para.mu: the parameter \mu in the paper
%       para.lamda: the parameter \gamma in the paper
% W: projection matrix
% F: soft label matrix

% Ref: Feiping Nie, Dong Xu, Ivor W. Tsang, Changshui Zhang.
%      Flexible Manifold Embedding: A Framework for Semi-supervised and Unsupervised Dimension Reduction}.
%      Accepted by IEEE Transactions on Image Processing (TIP)





[dim,n] = size(X);

Xm = mean(X,2);
Xc = X - Xm*ones(1,n);

if dim < n
    St = Xc*Xc';
    A = inv(St+1/para.lamda*eye(dim))*Xc;  
else
    K = Xc'*Xc;
    A = Xc*inv(K+1/para.lamda*eye(n));
end;
N = para.mu*para.lamda*Xc'*A;
Lc = eye(n) - 1/n*ones(n);

M = L + para.mu*para.lamda*Lc - N;
M = (M+M')/2;

[q,r] = qr(ones(n,1));
R = q(:,2:end);
Mr = R'*M*R;
Mr = (Mr+Mr')/2;

[ev, eigvalue] = eig(Mr);
eigval = diag(eigvalue);
[eval_b, iEvals] = sort(eigval);
F = R*ev(:,iEvals(1:dim));
W = A*F;

W = W./(ones(size(W,1),1)*sqrt(sum(W.*W)));









