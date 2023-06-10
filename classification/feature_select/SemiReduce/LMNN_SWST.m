function [ W ] = LMNN_SWST( XTrain, LTrain )
%LMNN_SWST Summary of this function goes here
%   Detailed explanation goes here
K = 3;
a = 0.1;
b= 0.1;

stdd = std(XTrain);
meandd = mean(XTrain);
stdd = repmat(stdd, size(XTrain, 1), 1);
meandd = repmat(meandd, size(XTrain, 1), 1);
XTrain = (XTrain - meandd)./stdd;

XTrain = XTrain';

[L,~] = lmnn2(XTrain, LTrain',2);
XTrain=L*XTrain;
distance = pdist2(XTrain', XTrain');
[sorted,index] = sort(distance);
neighborhood = index(:, 1:(1+K));
sita = sorted(K+1, :);
sita = sita' *  sita;
A = exp(-(distance.^2)./sita);
diagA = diag(ones(size(distance, 2), 1));
A = A - diagA;
M = zeros(size(A, 1), size(A, 2));
for ii = 1 : size(distance, 1)
    M(neighborhood(ii, :), ii) = A(neighborhood(ii, :), ii);
    M(ii, neighborhood(ii, :)) = A(ii, neighborhood(ii, :));
end
    
[Sb, Sw] = compute_swsb(XTrain, LTrain);

LV = Sb;% + b.* X * Lf * X';
RV = Sw + Sb + a.* XTrain * M * XTrain';
%[U1,D1,QP] = svd(inv(RV)*LV);
%[V,D]=eig(RV, LV);
[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end

end

function [Sb, Sw] = compute_swsb(X, labels)
X = X';
%fprintf(1,'compute_swsb running on %d points in %d dimensions\n',N,D);
[classes, bar, labels] = unique(labels);
nc = length(classes);
Sw = zeros(size(X, 2), size(X, 2));
    
% Compute total covariance matrix
St = cov(X);

	% Sum over classes
for i=1:nc
        
    % Get all instances with class i
    cur_X = X(labels == i,:);

	% Update within-class scatter
	C = cov(cur_X);
	p = size(cur_X, 1) / (length(labels) - 1);
	Sw = Sw + (p * C);
end
    
    % Compute between class scatter
Sb       = St - Sw;
Sb(isnan(Sb)) = 0; Sw(isnan(Sw)) = 0;
Sb(isinf(Sb)) = 0; Sw(isinf(Sw)) = 0;
end