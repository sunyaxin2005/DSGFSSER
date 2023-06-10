function [ W  ] = get_spare_graph( X,  K)
%SPARE_GRAPH Summary of this function goes here
%   Detailed explanation goes here
 %       X = pdist2(X, X);
%         [sorted,index] = sort(distance);
%         sita = sorted(K+1, :);
%         sita = sita' *  sita;
%         X = exp(-(distance.^2)./sita);
[N, D] = size(X);
dis = pdist2(X, X, 'cosine');
[sorted,index] = sort(dis);
neighborhood = index(:, 2:(1+K));
W = zeros(N, N);
tau = .35;
% set tolA
tolA = 1.e-3;
for ii = 1 : N
    y = X(ii, :)';
    A = X(neighborhood(ii, :), :)';
    [coef,coef_debias,obj_GPSR_Basic,times_GPSR_Basic,debias_s,mses_GPSR_Basic]= ...
        GPSR_Basic(y,A,tau,...
        'Debias',0,...
        'StopCriterion',1,...
        'ToleranceA',tolA);
    W(ii, neighborhood(ii, :)) = coef;
    W( neighborhood(ii, :), ii) = coef;
end
for ii = 1 : N
     W(ii, :) =   W(ii, :)./sum( W(ii, :));
end
% I = eye(N, N);
% W = (I- W)' * (I - W);
W = W'+ W - W' * W;
end
