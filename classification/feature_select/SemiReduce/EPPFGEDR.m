function [ W ] = EPPFGEDR(  XL, L, XUL )
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

%第一步：计算M
[M] = getLLEGRAPH(X,K);

[Sg, Sl] = GetWlWg( XL', L, X' );
Dg = sum(Sg(:,:),2);
Lg = diag(Dg)-Sg;
Dl = sum(Sl(:,:),2);
Ll = diag(Dl)-Sl;

[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
% [ sb, sw ] = get_sbsw( XL', L );
LB = XL * LB* XL';
Lg =  X * Lg * X';
LW = XL * LW* XL';
M =  X * M * X';
Ll = XL * Ll * XL';
E =  X * E * X';
LV = LB/trace(LB) + b.* Lg/trace(Lg);
RV =  Ll/trace(Ll) + b * E /trace(E );

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
K = 20;
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