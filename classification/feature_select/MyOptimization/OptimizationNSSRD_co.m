function [W] = OptimizationNSSRD_co(X,Y, Para)
%OPTIMIZATIONNSSRD 此处显示有关此函数的摘要
%   此处显示详细说明
[n, m] = size(X);
alpha = 0.1;
lamda = Para.lambda;
beta = Para.beta;
K = size(Y, 2);
distance = pdist2(X, X);
[sorted,index] = sort(distance, 2);
sita1 = sorted(:, floor(2 * K)+1);
sita = mean(mean(sita1));
Ws = exp(-distance.^2/sita.^2);

distance = pdist2(X', X');
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
XtX = X'* X;
XtCX = XtX+ alpha * Lp;
XtY = X'* Y;
Ut3 = eye(m, m);
Ut2 = eye(m, m);
Ut1 = eye(m, m);
Ut3U = eye(m, m);
I1 = eye(m, m);
I = eye(n, n);
UT = eye(m, m);
U2 = eye(m, m);

P2 = (Ut2*Ut3 * XtX+ lamda * Ut2 + beta * Ut3+ 0.00000001 * I1)^-1 * Ut2*Ut3 * XtY;
for t=1:15
   P = P .*(( X'*S+beta*Wp*P)./(XtX*P+beta*Dp*P+alpha*U*P + alpha*UT*P));
   S = S .*((X*P+beta *Ws * S+lamda*S)./(S+beta*Ds * S+lamda*S*S'*S));
   U = 1./(diag((sqrt(sum(P.^2, 2))).*2));
   
   P2 = P2 .*(( XtY)./(XtX*P2+alpha*U2*P2+ alpha*UT*P2));
   U2 = 1./(diag((sqrt(sum(P2.^2, 2))).*2));

%     P2 = (Ut1*Ut3 * XtX+ lamda * Ut1 + b*beta * Ut3 + 0.00000001 * I1)^-1 * Ut1*Ut3 * XtY1;
%     Ut1 = diag(sqrt(sum(W1.^2, 2)).*2);
   
   PT = [0.01*P, P2];
%    Ut3 = 
   UT = 1./(((sqrt(sum(PT.^2, 2))).*2));
   
   W = PT;
end
% for t=1:15
%     W = (Ut1*Ut3 * XtX+ lamda * Ut1 + beta * Ut3 + 0.00000001 * I1)^-1 * Ut1*Ut3 * XtY;
%     Ut1 = diag(sqrt(sum(W.^2, 2)).*2);
% %     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
% %     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
% end

% for t =  1 : 15
%     W1 = (Ut1*Ut3 * XtX+ 0.5*lamda * Ut1 + beta * Ut3 + 0.00000001 * I1)^-1 * Ut1*Ut3 * X'*S;
%     S = (4 * I + 2* Ls)^-1 * X*W1;
%     Ut1 = diag(sqrt(sum(W1.^2, 2)).*2);
%     
%     W2 = (Ut2*Ut3 * XtX+ lamda * Ut2 + beta * Ut3 + 0.00000001 * I1)^-1 * Ut2*Ut3 * XtY;
%     Ut2 = diag(sqrt(sum(W2.^2, 2)).*2);
%     
%     W = [0.5*W1, W2];    
%     Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
% end

end

