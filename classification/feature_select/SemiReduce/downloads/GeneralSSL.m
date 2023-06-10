% General Semi-supervised Learning with Graph
function F = GeneralSSL(A, T, alpha_l, alpha_u, is_normalize)
% A: n*n similarity matrix
% T: n*(c+1) class indicator matrix. Tij=1 if xi is labeled as j, Tij=0 otherwise
% alpha_l: parameter associated with labeled data points
% alpha_u: parameter associated with unlabeled data points
% is_normalize: if uses normalized weights
% F: n*(c+1) soft label matrix

% Ref: Changshui Zhang, Feiping Nie, Shiming Xiang. A General Kernelization
%      Framework for Learning Algorithms Based on Kernel PCA. Neurocomputing, 2010, 73(4-6): 959--967.




n = size(A,2);

labeled_idx = find(sum(T,2) == 1);
t = ones(n,1);
t(labeled_idx) = 0;
T = [T, t];


D = sum(A,2);
invD  = spdiags(1./D,0,n,n);
if is_normalize == 1 % normalized weights
    D1 = sqrt(invD);
    norA = D1*A*D1; 
    norD = sum(norA,2);
    P = spdiags(1./norD,0,n,n)*norA;
else
    P = invD*A;
end;

I_a = alpha_u*ones(n,1);I_a(labeled_idx) = alpha_l; I_a = spdiags(I_a,0,n,n);
I_b = speye(n) - I_a;
H = speye(n) - I_a * P;
F = H\(I_b*T);