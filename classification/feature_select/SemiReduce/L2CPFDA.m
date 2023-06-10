function [W] = L2CPFDA(XL, L, XUL)
%L2CPFDA 此处显示有关此函数的摘要
%   此处显示详细说明
%Dimensionality reduction by collaborative preserving Fisher discriminant
%analysis, neurocomputing, 
n = size(XL, 1);
[ S ] = get_spare_graph_L2( XL);
[WW, Wb] = compute_WW_Wb(XL, L);
I = eye(n, n);
LW = I - WW * S' - S* WW'+ S * S';
Lb = -Wb* S'-S*Wb';
LV = XL' * Lb * XL;
RV = XL' * LW * XL;


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


function [WW, Wb] = compute_WW_Wb(XL, L)
n = size(XL, 1);
UL = unique(L);
NC = zeros(length(UL), 1);
for ii = 1 : length(UL)
    NC(UL(ii)) = length(L(L==UL(ii)));
end
WW = zeros(n, n);
Wb = zeros(n, n);
for ii = 1 : n
    for jj =  1 : n
        if L(ii) == L(jj)
            WW(ii, jj) = 1/NC(L(ii));
            Wb(ii, jj) = 0;
        else
            WW(ii, jj) = 0;
            Wb(ii, jj) = 1/(n-NC(L(ii)));
        end
    end
end

end