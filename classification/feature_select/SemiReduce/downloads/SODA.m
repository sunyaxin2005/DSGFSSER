function W = SODA(X, T, A, para)
% X: n*dim matrix, each row is a data point
% T: n*c matrix, class indicator matrix. Tij=1 if xi is labeled as j, Tij=0 otherwise
% A: n*n matrix, similarity matrix
% para: regularization paramter and reduced dimension
% W: d*reduced_dim matrix, projection matrix

% Ref: Feiping Nie, Shiming Xiang, Yangqing Jia, Changshui Zhang.
%      Semi-Supervised Orthogonal Discriminant Analysis via Label Propagation. 
%      Pattern Recognition (PR), Volume 42, Issue 11, Pages 2615-2627, 2009. 




lambdascale = para.lambdascale;
reduced_dim = para.reduced_dim;

[n, dim] = size(X);

if isfield(para,'semiF') % if soft labels have been provided by other semi-supervised learning method
    F = semiF;
else
    F = GeneralSSL(A, T, 0, para.alpha_u, 0);
end;


FF = F(:,1:end-1);
B = diag(sum(FF,2));
D = diag(1./sum(FF,1));
a = sum(B(:));
Lt = B - 1/a*B*ones(n,1)*ones(1,n)*B;
Lw = B - FF*D*FF';

St = X'*Lt*X;
Sw = X'*Lw*X;
Sb = St - Sw;
lambda = lambdascale* max(diag(Sw));
Sw = Sw + lambda*eye(dim);

W = TraceRatio_fast(Sb, Sw, reduced_dim);




