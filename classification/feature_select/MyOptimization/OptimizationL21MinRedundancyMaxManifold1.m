function [ W  ] = OptimizationL21MinRedundancyMaxManifold1(  X, Y, M, L, Para)
%OPTIMIZATIONL21MANIFOLD Summary of this function goes here
%   Detailed explanation goes here
mean_x = mean(mean( X'* X));
mean_L = mean(mean(L));
miu = Para.beta * mean_x/mean_L;
lamda = Para.lambda;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:样本个数*样本维数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
I1 = eye(m, m);
Ut = eye(m, m);
XtY =  X'* Y;
XtCX = X'* (I + miu * M) * X ;
for t=1:1
    W = (Ut* XtCX + lamda * I1)^-1 * Ut * XtY;
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
end
H = zeros(m, m);
for t = 1 : 15
    for ii = 1 : m
        for jj = 1 : m
            H(ii, ii) = H(ii, ii) + norm(W(jj, :), 2) * L(ii, jj);
        end
        H(ii, ii) = H(ii, ii)/(m * (m - 1) * norm(W(jj, :), 2));
    end
    W = (Ut* XtCX + lamda * I1 + lamda * H)^-1 * Ut * XtY;
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
    
end
end

