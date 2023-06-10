function [ W  ] = ELDE(  XL, L, XUL )
%ELDE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%Matrix exponentialbasedsemi-superviseddiscriminantembedding for
%imageclassification Pattern recognition 2017
X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;

%��������
K = 1;
a = 0.1;
b = 0.1;
X = X';
XL = XL';
%��һ��������M
[M] = getLLEGRAPH(X,K);
%M = M.^2;
%�ڶ���������k��Զ������˹ͼ
[GF, GN] = get_graphnf(X, 6);
Df = sum(GF(:,:),2);
Lf = diag(Df)-GF;
Dn = sum(GN(:,:),2);
Ln = diag(Dn)-GN;
%�������������������ɢ������
[SB, SW, ST] = get_SWSB(XL, L);
DB = sum(SB(:,:),2);
LB = diag(DB)-SB;
DW = sum(SW(:,:),2);
LW = diag(DW)-SW;
% DT = sum(ST(:,:),2);
% LT = diag(DT)-ST;
% LV = XL * LB* XL' + b.* X * Lf * X';
% RV =  XL * LW* XL' + a.* X * M * X';
LB = XL * LB* XL';
Lf = X * Lf * X';
LW = XL * LW* XL';
Ln= X * Ln * X';
LV = LB/trace(LB);% + b.* Lf/trace(Lf)';
RV =  LW/trace(LW) + a.* Ln/trace(Ln);
%[U1,D1,QP] = svd(inv(RV)*LV);
[V,D]=eig(exp(LV), exp(RV));
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
V = diag(D);
for i = 1:size(W,2)
    W(:,i) = W(:,i)./norm(W(:,i));
end

end

