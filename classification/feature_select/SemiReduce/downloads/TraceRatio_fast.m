% Fast algorithm to solve the TraceRatio problems: 
% max_W   tr(W'*A*W)/tr(W'*B*W)
% min_W   tr(W'*A*W)/tr(W'*B*W)
% The matrix A should be symmetric and B must be positive semi-definite.
function [W, obj] = TraceRatio_fast(A, B, dim, isMax)
% Ref: Yangqing Jia, Feiping Nie, Changshui Zhang. Trace Ratio Problem Revisited. 
%      IEEE Transactions on Neural Networks (TNN), Volume 20, Issue 4, Pages 729-735, 2009.
%
% We give a fast algorithm in this paper to solve the TraceRatio problem.



if nargin < 4
    isMax = 1;
end;

n = size(A,1);
W = eye(n,dim);
obj = trace(W'*A*W)/trace(W'*B*W);

err = 1;
counter = 0;
while err > 10^-3 && counter < 20
    obj0 = obj;
    M = A - obj*B;  M = max(M,M');
    [v d] = eig(M);
    [W idx] = vector_selection(A, B, v, dim, isMax);

    obj = trace(W'*A*W)/abs(trace(W'*B*W));
    err = abs(obj - obj0)/abs(obj0);
    counter = counter + 1;
end;

if counter == 20
    disp('error! trace ratio can not converge');
end;



%% Select eigenvectors to maximize the trace ratio objective
function [W, idx] = vector_selection(A, B, v, dim, isMax)

n = size(A,2);
a = diag(v'*A*v);
b = abs(diag(v'*B*v));
idx = n-dim+1:n;
obj = sum(a(idx))/sum(b(idx));

err = 10;
counter = 0;
while err > 0 && counter < 20
    obj0 = obj;
    m = a - obj*b;
    if isMax == 1
        [temp idx] = sort(m, 'descend');
    else
        [temp idx] = sort(m);
    end;
    idx = idx(1:dim);

    obj = sum(a(idx))/sum(b(idx));
    err = abs(obj - obj0);
    counter = counter + 1;
end;

if counter == 20
    disp('error! trace ratio can not converge');
end;

W = v(:,idx);