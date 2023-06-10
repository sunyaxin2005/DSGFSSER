function [ W ] = SODA_main(  XL, L, XUL )
%SODA_MAIN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[GF, GN] = get_graphnf(XL', 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
[T] = get_matrix_label(L);

d = size(Dn, 1)/2;
if size(Dn, 1)/2 > 200
    d =200;
end
 para.lambdascale = 0.1;
 para.alpha_u = 0.1;
para.reduced_dim = d ;
W = SODA(XL, T, GN, para);
end

function [Y] = get_matrix_label(L)
c = length(unique(L));
Y = zeros(length(L), c);
for ii = 1 : length(L)
    Y(ii, L(ii)) = 1;
end
end