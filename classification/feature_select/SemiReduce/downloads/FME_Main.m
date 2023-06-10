function [  W ] = FME_Main(  XL, L, XUL  )
%FME_MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
[GF, GN] = get_graphnf(XL', 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
[T] = get_matrix_label(L);
 para.ul = size(XL, 1);
 para.uu = 0;
 para.mu = 0.1;
 para.lamda = 0.1;
 %       para.ul: parameter associated to labeled data points
%       para.uu: parameter associated to unlabeled data points(usually be 0)
%       para.mu: the parameter \mu in the paper
%       para.lamda: the parameter \gamma in the paper
[a, b, F] = FME_semi(XL',Ln, T, para);
W.a = a;
W.b = b;
end

function [Y] = get_matrix_label(L)
c = length(unique(L));
Y = zeros(length(L), c);
for ii = 1 : length(L)
    Y(ii, L(ii)) = 1;
end
end
