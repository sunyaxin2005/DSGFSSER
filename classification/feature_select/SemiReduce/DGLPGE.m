function [W] = DGLPGE(XL, L, XUL)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%Discriminative globality and locality preserving graph emb e dding for dimensionality reduction
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

b = 0.1;
[W_plus, W_jian, W_n] = compute_W_plus(XL, L);
DB = sum(W_jian(:,:),2);
LB = diag(DB)-W_jian;
DW = sum(W_plus(:,:),2);
LW = diag(DW)-W_plus;

DT = sum(W_n(:,:),2);
LT = diag(DT)-W_n;

LV = XL'*LB*XL ;
RV =  XL'*LW*XL + b * XL'*LT* XL;

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
function [W_plus, W_jian, W_n] = compute_W_plus(XL, L)

n = size(XL, 1);
m = size(XL, 2);

dis = pdist2(XL, XL);
UL = unique(L);

% t_plus = zeros(length(UL), 1);
% for ii = 1 : length(UL)
%     t_plus(ii) = 
% end

W_plus = zeros(n, n);
for ii = 1 : n
        for jj = 1 : n
        if(L(ii) == L(jj))
            aaa = exp(-dis(ii, jj).^2/sum(dis(ii, L==UL(L(ii))).^2)* sum(L==UL(L(ii))));
            W_plus(ii, jj) = aaa * (1 + aaa)/2;
        else
            W_plus(ii, jj) = 0;
        end
    end
end
W_plus = (W_plus+ W_plus')/2;

W_jian = zeros(n, n);
for ii = 1 : n
    for jj = 1 : n
        if(L(ii) ~= L(jj))
            aaa = exp(-dis(ii, jj).^2/sum(dis(ii, L~=UL(L(ii))).^2) * sum(L~=UL(L(ii))));
            W_jian(ii, jj) = aaa * (1 + aaa)/2;
        else
            W_jian(ii, jj) = 0;
        end
    end
end
W_jian = (W_jian+ W_jian')/2;

[sorted,index] = sort(dis);
K = 5;
for i=1:n
    for j=2:1+K
        aaa = exp(-dis(i, index(i,j)).^2/sum(dis(ii, index(i,2:K))) * K);
        W_n(i,index(i,j)) = aaa * (1 + aaa)/2;
        W_n(index(i,j),i) = aaa * (1 + aaa)/2;
    end
end

end
