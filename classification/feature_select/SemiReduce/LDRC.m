function [ W ] = LDRC( XL, L, XUL  )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
ata = 0.1;
XL = XL';
XUL = XUL';
C = unique(L);
C = length(C);
H = cell(C, 1);
for ii = 1 : C
    Xii = XL(:, L == ii);
    H{ii} = Xii * (Xii'*Xii)^-1 * Xii';
end

Xinter = cell(C, C);
Xintra = cell(C, 1);
for ii = 1 : C
    Xii = XL(:, L == ii);
    for jj = 1 : C
        Xinter{ii, jj} = H{jj} * Xii;
    end
    Xintra{ii} = H{ii} * Xii;
end
Eb = 0;
for ii = 1 : C
    Xii = XL(:, L == ii);
    for jj = 1 : C
        if ii == jj
            continue;
        end
        Eb = Eb+(Xii - Xinter{ii, jj}) * (Xii - Xinter{ii, jj})';
    end
end
Eb = Eb/(size(XL, 2) * (C-1));
Ew = 0;
for ii = 1 : C
    Xii = XL(:, L == ii);
    Ew = (Xii - Xintra{ii}) * (Xii - Xintra{ii})';
end
Ew = Ew/size(XL, 2) ;
I = eye(size(Ew));
LV = Eb;
RV = Ew + ata * I;
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

