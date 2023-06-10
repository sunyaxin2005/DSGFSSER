% max_{F,W'*W=I} trace(F'*Lb*F) / (trace(F'*La*F)+||X'*W-F||^2+r*||W||^2)
% La should be positive semi-definite
% La-Lb should be positive semi-definite
function [W, F, lambda, converge] = TR_FSDA(X, Lb, La, r, m)
% X: d*n training data matrix
% Lb, La: n*n Laplacian matrices defined as in the paper
% r: regularization parameter
% m: reduced dimensionality
% W: d*m projection matrix
% F: n*m matrix
% lambda: objective values in the iterations
% converge: the value should near 0 iff the algorithm converges

% Ref: Yi Huang, Dong Xu, Feiping Nie.
%      Semi-supervised Dimension Reduction using Trace Ratio Criterion.
%      IEEE Transactions on Neural Networks and Learning Systems, 23(3):519-526, 2012. 

[d n] = size(X);

W = orth(rand(d,m));
F = X'*W;
lam0 = trace(F'*Lb*F) / (trace(F'*La*F)+r*trace(W'*W));
lam1 = 1;
lambda(1) = lam0;
for iter = 1:20
    lam = (lam0+lam1)/2;
    M = lam*La-Lb + lam*eye(n);
    [v dd] = eig(M); dd = diag(dd); mind = min(dd);
    if mind <= 0
        lam0 = lam;
    else
        M = X*(eye(n)-lam*inv(lam*La-Lb+lam*eye(n)))*X'+r*eye(d);
        [v dd] = eig1(M,m,0); obj = sum(dd);
        if obj > 0
            lam1 = lam;
        else
            break;
        end;
    end;
end;
iter

lambda(2) = lam;
for iter = 1:10
    lam = lambda(iter+1);
    inL = inv(lam*La-Lb + lam*eye(n));
    M = X*(eye(n)-lam*inL)*X'+r*eye(d);
    [W dd] = eig1(M,m,0);
    F = lam*inL*X'*W;
    E = X'*W-F;
    lambda(iter+2) = trace(F'*Lb*F) / (trace(F'*La*F)+trace(E'*E)+r*trace(W'*W));
end;
converge = sum(dd);






function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)

if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end;

if nargin < 3
    isMax = 1;
    isSym = 1;
end;

if nargin < 4
    isSym = 1;
end;

if isSym == 1
    A = max(A,A');
end;
[v d] = eig(A);
d = diag(d);
%d = real(d);
if isMax == 0
    [d1, idx] = sort(d);
else
    [d1, idx] = sort(d,'descend');
end;

idx1 = idx(1:c);
eigval = d(idx1);
eigvec = v(:,idx1);

eigval_full = d(idx);