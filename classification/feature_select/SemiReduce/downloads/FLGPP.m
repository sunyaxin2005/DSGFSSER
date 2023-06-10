% min_{F,W'*W=I} (trace(F'*L*F) + r*||X'*W-F||^2) / trace(F'*Lp*F)
% Lp should be positive definite
function [W, F, converge] = FLGPP(X, L, Lp, r, m)
% X: d*n training data matrix, each column is a data point
% L, Lp: two Laplacian matrices
% r: parameter
% m: reduced dimension
% W: d*m projection matrix
% F: n*m embedding matrix
% converge: if converge -> 0, then the algorithm converges

% Ref:
% Feiping Nie, Xiao Cai, Heng Huang. 
% Flexible Shift-Invariant Locality and Globality Preserving Projections. 
% The European Conference on Machine Learning and Principles and Practice of Knowledge Discovery in Databases (ECML PKDD), Nancy, France, 2014.


[d n] = size(X);

W = orth(rand(d,m));
F = X'*W;
lam1 = trace(F'*L*F) / trace(F'*Lp*F);
lam0 = 0;
lambda(1) = lam1;
for iter = 1:20
    lam = (lam0+lam1)/2;
    M = L-lam*Lp + r*eye(n);
    [v dd] = eig(M); dd = diag(dd); mind = min(dd);
    if mind <= 0
        lam1 = lam;
    else
        M = X*(eye(n)-r*inv(L-lam*Lp+r*eye(n)))*X';
       %[v dd] = eig1(M,m,0); obj = sum(dd);
        [V,D] = mySVD(M);
        [sorted,index] = sort(diag(D),'descend');
        V1=V;
        for i=1:size(index)
            V1(:,i) = V(:,index(i));
        end
        if m > size(D, 1)
            m = size(D, 1);
        end
        dd = D(1:m, 1:m);
        obj = sum(dd);
        if obj > 0
            lam0 = lam;
        else
            break;
        end;
    end;
end;
iter

lambda(2) = lam;
for iter = 1:10
    inL = inv(L-lambda(iter+1)*Lp+r*eye(n));
    M = X*(eye(n)-r*inL)*X';
    %[W dd] = eig(M);%eig1(M,m,0);
            [V,D] = mySVD(M);
        [sorted,index] = sort(diag(D),'descend');
        V1=V;
        for i=1:size(index)
            V1(:,i) = V(:,index(i));
        end
                if m > size(V1, 2)
            m = size(V1, 2);
        end
        W = V1(:, 1:m);
        dd = sorted(1:m);
    lambda(iter+2) = sum(dd) / (r*trace(W'*X*inL*Lp*inL*X'*W))+lambda(iter+1);
end;
F = r*inL*X'*W;
converge = sum(dd);