function [ W ] = LRGAAFFSRC(  XL, Label, XUL, K)
%LRGAAFF Summary of this function goes here
%   Detailed explanation goes here
X = [XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;
X = X';
XL = XL';
N = size(X, 2);
%定义U，按离第rank_index个元素的远近进行排序，那么第rank_index对角的值赋值非常大，其它为1
lambuda = 1000;
a = 0.1;
b = 0.1;
gamma = 0.5;
distance = pdist2(X', X');
[sorted,index] = sort(distance);
neighborhood = index(:, 1:(1+K));
H = zeros(K + 1, 1);
H = H + 1;
I = eye(K+1);
H = I - H * H';

S = [];
for ii = 1 : N
    Si = zeros(N, K + 1); 
    for jj = 1 : K + 1
        Si(neighborhood(ii, jj), jj) = 1;
    end
    S = [S, Si];
end

A = zeros(N * (K + 1),  N * (K + 1));
for ii = 1 : N
    Li = lambuda * H *(H * X(:, neighborhood(ii, :))' * X(:, neighborhood(ii, :)) * H + lambuda * I)^(-1) * H;
    start_i = (ii - 1) * (K + 1) + 1;
    end_i = start_i + (K + 1) - 1;
    A(start_i:end_i,  start_i:end_i) = Li;
end

L = S * A * S';
for ii = 1 : N
    L(ii, :) = L(ii, :)./sum( L(ii, :));
end
[SB, SW, ST] = compute_swsb(XL, Label);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
LV = XL * LB* XL';% + b.* X * Lf * X';
RV =  XL * LT* XL' + a.* X * L * X';
[V,D]=eig(LV, RV);

[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
end

function [SB, SW, ST] = compute_swsb(XL, L)
K = 5;
X = XL';
[N, D] = size(X);
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
A = exp(-(distance.^2)./sita);
diagA = diag(ones(size(distance, 2), 1));
A = A - diagA;
neighborhood = index(2:(1+K),:);
A1 = A;
for ii = 1 : size(neighborhood, 1)
    A1(ii, neighborhood(ii, :)) = 0; 
    A1(neighborhood(ii, :), ii) = 0; 
end
A1 = A - A1;
SW = zeros(N, N);
UL = unique(L);
nClass = length(UL);
A = ones(N,N);
%赋值有label的
for ii = 1 : nClass
    SW(L == UL(ii), L == UL(ii)) = A(L == UL(ii), L == UL(ii))/length(L(L == UL(ii)));
end% %赋值没有label的
%赋值没有label的
for ii = 1 : N
    SW(neighborhood(:, ii), ii) = SW(neighborhood(:, ii), ii) + A(neighborhood(:, ii), ii)/N;
    SW(ii, neighborhood(:, ii)) = SW(ii, neighborhood(:, ii)) + A(ii, neighborhood(:, ii))/N;
end
SW=SW-diag(diag(SW));
%赋值有label的
SB = zeros(N, N);
for ii = 1 : nClass
    SB(L == UL(ii), L == UL(ii)) = -A(L == UL(ii), L == UL(ii))*(1/length(L(L == UL(ii))) - 1/N);
    SB(L == UL(ii), L ~= UL(ii)) = A(L == UL(ii), L ~= UL(ii))/N;
    SB(L ~= UL(ii), L == UL(ii)) = A(L ~= UL(ii), L == UL(ii))/N;
end% %赋值没有label的
SB=SB-diag(diag(SB));
ST = A;
end

% function [Sb, Sw] = compute_swsb(X, labels)
% [D,N] = size(X);
% %fprintf(1,'compute_swsb running on %d points in %d dimensions\n',N,D);
% X = X';
% [classes, bar, labels] = unique(labels);
% nc = length(classes);
% Sw = zeros(size(X, 2), size(X, 2));
%     
% % Compute total covariance matrix
% St = cov(X);
% 
% 	% Sum over classes
% for i=1:nc
%         
%     % Get all instances with class i
%     cur_X = X(labels == i,:);
% 
% 	% Update within-class scatter
% 	C = cov(cur_X);
% 	p = size(cur_X, 1) / (length(labels) - 1);
% 	Sw = Sw + (p * C);
% end
%     
%     % Compute between class scatter
% Sb       = St - Sw;
% Sb(isnan(Sb)) = 0; Sw(isnan(Sw)) = 0;
% Sb(isinf(Sb)) = 0; Sw(isinf(Sw)) = 0;
% end
