function [ S ] = get_DSPNE(  XL, L )
%GET_DSPNE Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(XL);
% tau = 0.35;
% % set tolA
% tolA = 1e-5;
param.lambda=0.1;
param.L = 2;
param.numThreads=-1; % number of threads
XL = XL';
XL =XL ./repmat(sqrt(sum(XL .^2)),[size(XL ,1) 1]);
S = zeros(N, N);
for ii = 1 : size(XL, 2)
    inx = L == L(ii);
    inx(ii) = 0;
    A = XL(:, inx);
    y = XL(:, ii);
    
%     [coef,coef_debias,obj_GPSR_Basic,times_GPSR_Basic,debias_s,mses_GPSR_Basic]= ...
%         GPSR_Basic(y, A,tau,...
%         'Debias',0,...
%         'StopCriterion',1,...
%         'ToleranceA',tolA);
%    coefs1=mexLasso(y,A,param);
    [coefs, funVal]= nnLeastR(A, y, param.lambda);
    S(ii, inx) = coefs;
%    S(ii, ii) = 1; 
end
%S = S + S'/2;
% for ii = 1 : N
%      S(ii, :) =   S(ii, :)./sum( S(ii, :));
% end
I = eye(N, N);
S = I - S - S' + S'*S;
end

