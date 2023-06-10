function [ W ] = SWLDRC( XL, L, XUL  )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X = [XL;XUL];
%X = XL;
K = 5;
ata = 0.1;
XL = XL';
X = X';
C = unique(L);
C = length(C);
H = cell(C, 1);
for ii = 1 : C
    Xii = XL(:, L == ii);
    H{ii} = Xii * (Xii'*Xii)^-1 * Xii';
end


Eb = 0;
for ii = 1 : size(XL, 2)
    xi = XL(:, ii);
    li = L(ii);
    for jj = 1 : C
        if jj == li
            continue;
        end
        xijinter = H{jj} * xi;
        Eb = Eb + (xi - xijinter) * (xi - xijinter)' * exp(-pdist2(xi', xijinter'));
    end
end
Eb = Eb/(size(XL, 2) * (C-1));
Ew = 0;
for ii = 1 : size(XL, 2)
    xi = XL(:, ii);
    li = L(ii);
    xiintra = H{li} * xi;
    Ew = Ew + (xi - xiintra) * (xi - xiintra)' * exp(-pdist2(xi', xiintra'));
end
Ew = Ew/size(XL, 2) ;

X2 = sum(X.^2,1);
N = size(X, 2);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*(X'*X);
[sorted,index] = sort(distance);
Em = 0;
for ii = 1 : size(X, 2)
    XNii = X(:, index(2:K+1, ii));
    HXii = XNii * (XNii'*XNii)^-1 * XNii';
    xi = X(:, ii);
    xiintra = HXii * xi;
    Em = Em + (xi - xiintra) * (xi - xiintra)' * exp(-pdist2(xi', xiintra'));
end
Em = Em/(size(X, 2) * K) ;

I = eye(size(Ew));
LV = Eb/trace(Eb);
RV = (Ew + ata * I)/trace(Ew + ata * I) + ata * Em/trace(ata * Em);
[V,D]=eig(LV, RV);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end

% for i = 1:size(W,2)
%     W(:,i) = W(:,i)./norm(W(:,i));
% end
end

