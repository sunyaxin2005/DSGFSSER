function [ W ] = DOEPP(  XL, L, XUL )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X =[XL;XUL];
N = size(X, 1);
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;
K = 1;
a = 0.1;
b = 0.1;
r = 0.1;
u = 0.1;
X = X';
XL = XL';



yi = zeros(N, 1);
I = eye(N);
Hc = I - 1/N * (yi*yi');
Xc = X*Hc;
A = r*(r*(Xc*Xc')+eye(size(Xc*Xc')))^-1*Xc;
B = Hc*X'*A+1/N * (yi*yi');
E = u*(A'*A)+u*r*(B-I)'*(B-I);

[Dl, Ml] = GetMD(XL, L);
[SB, SW, ST] = get_SWSB(XL, L);
DW = sum(ST(:,:),2);
L = diag(DW)-ST;
L1 = L + r * Ml;

Dl = XL * Dl* XL';
L1 = XL * L1* XL';
LV = Dl/trace(Dl);
RV = L1/trace(L1);

% LV = XL * LB* XL' + b.* X * Lg * X';
% RV =  XL * LW* XL' + a.* X * M * X';

% LV = sb + b.* X * Lg * X';
% RV =  sw + a.* XL * Ll * XL';

%[U1,D1,QP] = svd(inv(RV)*LV);
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

function [Wg, Wl] = GetWlWg( XL, L, X )
UL = unique(L);
K = 10000;
for jj = 1 : length(UL)
   if ( K > length(L(L == jj)))
       K = length(L(L == jj));
   end
end
K = 2 * K;
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
t = mean(mean(distance));
Wg = distance.^2 * exp(-distance.^2./(2 * t.^2));
distance = pdist2(XL, XL);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
t = mean(mean(distance));
Wl = exp(-distance.^2./(2 * t.^2));
for ii = 1 : size(Wl, 1)
    for jj = 1 : size(Wl, 2)
        if(L(ii) ~= L(jj))
            Wl(ii, jj) = 0;
        end
    end
end
end

function [Dl, Ml] = GetMD(XL, L)
K = 5;
X = XL';
[N, D] = size(X);

% distance = pdist2(X, X);
% [sorted,index] = sort(distance);
% sita = sorted(K+1, :);
% sita = sita' *  sita;
% A = exp(-(distance.^2)./sita);
% diagA = diag(ones(size(distance, 2), 1));
% A = A - diagA;
% neighborhood = index(2:(1+K),:);
A = ones(N,N);
SW = zeros(N, N);
UL = unique(L);
nClass = length(UL);
%赋值有label的
% for ii = 1 : nClass
%     SW(L == UL(ii), L == UL(ii)) = A(L == UL(ii), L == UL(ii))/length(L(L == UL(ii)));
% end% %赋值没有label的
for ii = 1 : N
    for jj = 1 : N
        if(L(ii) == L(jj)) && ii ~= jj
            SW(ii, jj) = 1/length(L(L==L(jj)));
        end
    end
    %SW(ii, :) = SW(ii, :)./sum(SW(ii, :));
end
%赋值没有label的
% for ii = 1 : N
%     SW(neighborhood(:, ii), ii) = SW(neighborhood(:, ii), ii) + A(neighborhood(:, ii), ii)/N;
%     SW(ii, neighborhood(:, ii)) = SW(ii, neighborhood(:, ii)) + A(ii, neighborhood(:, ii))/N;
% end
SW=SW-diag(diag(SW));
%赋值有label的
SB = zeros(N, N);
% for ii = 1 : nClass
%     SB(L == UL(ii), L == UL(ii)) = -A(L == UL(ii), L == UL(ii))*(1/length(L(L == UL(ii))) - 1/N);
%     SB(L == UL(ii), L ~= UL(ii)) = A(L == UL(ii), L ~= UL(ii))/N;
%     SB(L ~= UL(ii), L == UL(ii)) = A(L ~= UL(ii), L == UL(ii))/N;
% end% %赋值没有label的
for ii = 1 : N
    for jj = 1 : N
        if(L(ii) ~= L(jj)) && ii ~= jj
            
            SB(ii, jj) = 1/(N - length(L(L==L(jj))));
        end
    end
     %SB(ii, :) =  SB(ii, :) ./sum(SB(ii, :));
end
I = eye(size(SB));
SB=SB-diag(diag(SB));
DB = diag(sum(SB));
DW = diag(sum(SW));
Ml = 3 * I + DB + SB + SB' - 2*SW;
Dl = 1 + DB;
end