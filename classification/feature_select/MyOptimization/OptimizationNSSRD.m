function [P] = OptimizationNSSRD(X,Y, Para)
%OPTIMIZATIONNSSRD 此处显示有关此函数的摘要
%   此处显示详细说明
m = size(X, 2);
X = X';
alpha = 0.1;
lamda = Para.lambda;
beta = Para.beta;
K = size(Y, 2);
distance = pdist2(X', X');
[sorted,index] = sort(distance, 2);
sita1 = sorted(:, floor(2 * K)+1);
sita = mean(mean(sita1));
Ws = exp(-distance.^2/sita.^2);

distance = pdist2(X, X);
[sorted,index] = sort(distance, 2);
sita1 = sorted(:, floor(2 * K)+1);
sita = mean(mean(sita1));
Wp = exp(-distance.^2/sita.^2);

Ds = diag(sum(Ws(:,:),2));
Ls = (Ds)-Ws;
Dp = diag(sum(Wp(:,:),2));
Lp = (Dp)-Wp;


S = Y;
U = eye(m, m);
[V,D] = eig(Lp);
[sorted,index] = sort(diag(D),'descend');
for i=1:size(S, 2)
    P(:,i) = V(:,index(i));
end
XtX = X* X' ;
for t=1:15
   P = P .*(( X*S+beta*Wp*P)./(XtX*P+beta*Dp*P+alpha*U*P));
   S = S .*((X'*P+beta *Ws * S+lamda*S)./(S+beta*Ds * S+lamda*S*S'*S));
   U = 1./((diag(sqrt(sum(P.^2, 2)).*2)));
end


end

