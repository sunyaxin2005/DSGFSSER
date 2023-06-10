function [ W ] = FLGPP_main( XL, L, XUL  )
%FLGPP_MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
% [SB, SW, ST] = get_SWSB(XL', L);
% DB = sum(SB(:,:),2);
% LB = diag(DB)-SB;
% DW = sum(SW(:,:),2);
% LW = diag(DW)-SW;
% DW = diag(DW);
% I1 = zeros(size(DW, 1), 1) + 1;
% Lp = DW - DW * (I1 * I1') * DW/(I1'*DW*I1);
[GF, GN] = get_graphnf(XL', 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
 Dn = diag(Dn);
I1 = zeros(size(Dn, 1), 1) + 1;
Lp = Dn - Dn * (I1 * I1') * Dn/(I1'*Dn*I1);

d = size(Dn, 1)/2;
if size(Dn, 1)/2 > 200
    d =200;
end
[W, F, converge] = FLGPP(XL', Ln, Lp, 0.1, d);

end

