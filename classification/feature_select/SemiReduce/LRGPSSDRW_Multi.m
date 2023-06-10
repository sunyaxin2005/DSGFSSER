function [ Ws ] = LRGPSSDRW_Multi(  XL, L, XUL, num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
     num = 3;
end

X =[XL;XUL];
stdd = std(X);
meandd = mean(X);
stdd1 = repmat(stdd, size(X, 1), 1);
meandd1 = repmat(meandd, size(X, 1), 1);
X = (X - meandd1)./stdd1;
stdd = repmat(stdd, size(XL, 1), 1);
meandd = repmat(meandd, size(XL, 1), 1);
XL = (XL - meandd)./1;


%参数设置
    K = 5;
    a = 0.1;
    b = 0.1;
    X = X';
    XL = XL';


    %第一步：计算M
    %[M] = getLLEGRAPH(X,K);
    [M] = get_NewLLEG(X',K);
    %M = M.^2;
    %第二步：计算k近远拉普拉斯图
    [GF, GN] = get_graphnf(X, 5);
    Df = sum(GF(:,:),2);
    Lf = diag(Df)-GF;
    Dn = sum(GN(:,:),2);
    Ln = diag(Dn)-GN;
    %第三步：计算类内类间散布矩阵
    [SB, SW, ST] = get_SWSB(XL, L);
    [WSB, WSW, WST] = get_weight_SW_SB(XL, L);
    %WSB = (WSB + SB)/2;
    %WSW = (WSW + SW)/2;
    DB = sum(SB(:,:),2);
    LB = diag(DB)-SB;
    DW = sum(SW(:,:),2);
    LW = diag(DW)-SW;
    DT = sum(ST(:,:),2);
    LT = diag(DT)-ST;

    WDB = sum(WSB(:,:),2);
    WLB = diag(WDB)-WSB;
    WDW = sum(WSW(:,:),2);
    WLW = diag(WDW)-WSW;
    
    LV = b.* X * Lf * X' +   XL * WLB * XL';
     RV =  a.* X * M * X' +   XL * WLW * XL';
%num = 3;
Ws = cell(num, 1);    
WiWiT = 0;
aa =  trace(XL * WLW * XL');
for ii = 1 :num    
   a1 = 0;
   if WiWiT ~= 0
       a1 = aa /trace(WiWiT) * 1; 
   end
   RV1 = RV + WiWiT * a1;
    %[U1,D1,QP] = svd(inv(RV)*LV);
    [V,D]=eig(LV, RV1);
    %[V,D,QP] = svd(pinv(RV1)*LV);
    [sorted,index] = sort(diag(D),'descend');
    W=V;
    for i=1:size(index)
        W(:,i) = V(:,index(i));
    end
     WiWiT = WiWiT + W*W';
    for i = 1:size(W,2)
        W(:,i) = W(:,i)./norm(W(:,i));
    end
    Ws{ii} = W;
end

end

