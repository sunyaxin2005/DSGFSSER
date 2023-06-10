function [ Ws ] = LRGPSSDRW_Multi_subspace(  XL, L, XUL, num )
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

    Ms = cell(num,1);
    numFeatures = size(X, 1);
    ratio = 1;

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
    
    Ws = cell(num, 1);    
    for ii = 1 : num
        usedFeatures = randperm(numFeatures, floor(numFeatures*ratio));
        %第一步：计算M
        %[M] = getLLEGRAPH(X,K);
        X_sub = X(usedFeatures, :);
        X_sub_L = XL(usedFeatures, :);
        [M] = get_NewLLEG(X_sub',K);
        Ms{ii} = M;
        [WSB, WSW, WST] = get_weight_SW_SB(X_sub_L, L);
            WDB = sum(WSB(:,:),2);
        WLB = diag(WDB)-WSB;
    WDW = sum(WSW(:,:),2);
    WLW = diag(WDW)-WSW;
        [GF, GN] = get_graphnf(X_sub, 5);
    Df = sum(GF(:,:),2);
    Lf = diag(Df)-GF;
    Dn = sum(GN(:,:),2);
    Ln = diag(Dn)-GN;
        
        LV = b.* X * Lf * X' +   XL * WLB * XL';
        RV =  a.* X * M * X' +   XL * WLW * XL';
        [V,D]=eig(LV, RV);
        W=V;
        [sorted,index] = sort(diag(D),'descend');
        for i=1:size(index)
            W(:,i) = V(:,index(i));
        end
        for i = 1:size(W,2)
             W(:,i) = W(:,i)./norm(W(:,i));
        end
        Ws{ii} = W;
    end
end

