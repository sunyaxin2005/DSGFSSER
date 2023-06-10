function [ W,C ] = SBDA_main(  XL, L, XUL )
%SBDA_MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
if(size(XL, 2) == 1024)
    [XL] = vector2symvector( XL );
    [XL] = vector2symvector( XL );
end
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
DT = sum(ST(:,:),2);
LT = diag(DT)-ST;
% [ sb, sw ] = get_sbsw( XL', L );
LB = XL * LB* XL';
Lg =  X * Lg * X';
LW = XL * LW* XL';
M =  X * M * X';
Ll = XL * Ll * XL';
E =  X * E * X';
LT = XL * LT * XL';
options.miu = 0.001;
options.card = 200;
options.lambda = 0.001;
options.LagMutiplier = 0.001;
options.ite = 50;
[W,C]=SBDA(LB,LW,LT,options);

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
