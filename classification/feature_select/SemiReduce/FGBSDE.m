function [ W ] = FGBSDE(  XL, L, XUL )
%FG 此处显示有关此函数的摘要
%   此处显示详细说明
r = 0.1;
u = 0.1;
a = 1;
b = 0.1;
lamda = 0.1;
X = [XL;XUL];
N = size(X, 1);
NL = size(XL, 1);

X = X';
XL = XL';

yi = zeros(N, 1);
I = eye(N);
Hc = I - 1/N * (yi*yi');
Xc = X*Hc;
A = r*(r*(Xc*Xc')+eye(size(Xc*Xc')))^-1*Xc;
B = Hc*X'*A+1/N * (yi*yi');

E = u*(A'*A)+u*r*(B-I)'*(B-I);
[ SB, SW, ST ] = get_SWSB( XL, L );
[GF, GN] = get_graphnf(X, 6);
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;

DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;

ML = 3 * eye(NL) + diag(DB) + SB +SB' - 2 * SW;
ML_ = zeros(N);
ML_(1:NL, 1: NL) = ML ;

DL = eye(length(DB)) + diag(DB);
DL_ = zeros(N);
DL_(1:NL, 1: NL) = DL ;

L1 = Ln + lamda * ML_;

[WSB, WSW, WST] = get_weight_SW_SB(XL, L);
WDB = sum(WSB(:,:),2);
WLB = diag(WDB)-WSB;
WDW = sum(WSW(:,:),2);
WLW = diag(WDW)-WSW;
% 
% LV = L1 + E;
% RV = DL_;
LV = XL * WLB* XL' + b.* X * Lf * X';
RV =  XL * WLW* XL' + a.* X * E * X';
[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end
% Z = W;
% W = A * Z;
% b = (Z'* yi - W'* X * yi);
% W = [W;b'];
end

