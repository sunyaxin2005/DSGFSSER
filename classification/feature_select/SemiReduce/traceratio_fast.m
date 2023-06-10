% objective: trace(W'*A*W)/trace(W'*B*W), s.t. W'*W = I;
function [W, lambda] = traceratio_fast(A, B, ReducedDim, isMax)
% isMax: 1, max trace(W'*A*W)/trace(W'*B*W)
%        0, min trace(W'*A*W)/trace(W'*B*W)

if nargin < 4
    isMax = 1;
end;
d = size(A,1);
W = eye(d,ReducedDim);
for iter = 1:10
    lambda(iter) = trace(W'*A*W)/trace(W'*B*W);
    M = A - lambda(iter)*B;
    M = (M+M')/2;
    [eigVec eigVal] = eig(M); eigVal = diag(eigVal);
    idx = FeatureSelect(A, B, eigVec, eigVal, ReducedDim, isMax);
    W = eigVec(:,idx(1:ReducedDim));
end;


function idx = FeatureSelect(A, B, W, eigVal, ReducedDim, isMax)

d = size(A,1);

a = diag(W'*A*W);
b = diag(W'*B*W);

if isMax == 1
    [ev, idx] = sort(eigVal,'descend');
else
    [ev, idx] = sort(eigVal,'descend');
end;
%idx = randperm(d); idx = idx';

idx_old = idx; idx_old(1:2) = idx_old(2:-1:1);
while sum(abs(idx - idx_old)) ~= 0
    idx_old = idx;
    lambda = sum(a(idx(1:ReducedDim)))/sum(b(idx(1:ReducedDim)));
    c = a - lambda*b;
    if isMax == 1
        [temp, idx] = sort(c,'descend');
    else
        [temp, idx] = sort(c);
    end;
end;
    




