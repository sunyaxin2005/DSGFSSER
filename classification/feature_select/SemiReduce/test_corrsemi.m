function [ prs ] = test_corrsemi(  X, L, ALg_type, sel_dims )
%TEST_CORRSEMI Summary of this function goes here
%   Detailed explanation goes here
ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'MMP', 'locateLRGPSSDR', 'LRGAAFF'};
random_times = 10;
X = NormalizeFea(X);
stdd = std(X);
meandd = mean(X);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;

load('E:\manifod_learning\caideng\test_data\新建文件夹\YaleB_32x32.mat');
gnd = NormalizeFea(fea);
stdd = std(gnd);
meandd = mean(gnd);
stdd = repmat(stdd, size(gnd, 1), 1);
meandd = repmat(meandd, size(gnd, 1), 1);
gnd = (gnd - meandd)./stdd;

pknn = zeros(random_times, length(sel_dims));
ptrasrc = zeros(random_times, length(sel_dims));
pkccsrc = zeros(random_times, length(sel_dims));
pkccsrc1 = zeros(random_times, length(sel_dims));
for ii = 1 : random_times
    try
    l = 6;
    ul = 6;
    us = 7;
    [train, semitrain, test] = get_l_ul_us_inx(L, l, ul);
    XLTrain = X(train, :);
    LTrain = L(train);
    XULTrain = X(semitrain, :);
    XTest = X(test, :);
    LTest = L(test);
 

    [eigvector, W] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, ALg_type);
    for jj = 1 : length(sel_dims)
        sel_dim = sel_dims(jj);
        if sel_dim > size(W, 2)
            sel_dim = size(W, 2);
        end
        XLTrain1 = XLTrain * eigvector * W(:, 1 : sel_dim);
        XULTrain1 = XULTrain * eigvector * W(:, 1 : sel_dim);
        XTest1 = XTest * eigvector * W(:, 1 : sel_dim);  
    
        [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
        p = length(LTest(classes == LTest))/length(classes);
        pknn(ii, jj) = p;
%         [ classes ] = SRC_MAIN( XLTrain, LTrain, XTest, 'CLSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc(ii, jj) = p;
        [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_src');
        p = length(LTest(classes == LTest))/length(classes);
        ptrasrc(ii, jj) = p;
        [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'kcckneighbor_src');
        p = length(LTest(classes == LTest))/length(classes);
        pkccsrc1(ii, jj) = p;
        
        if sel_dims(jj) > size(W, 2)
            break; 
        end
    end
    catch
    end
end
pknn = mean(pknn, 1);
pkccsrc = mean(pkccsrc, 1);
ptrasrc = mean(ptrasrc, 1);
pkccsrc1 =  mean(pkccsrc1, 1);
prs = [pknn; pkccsrc; ptrasrc; pkccsrc1];
end

function [train, semitrain, test] = get_l_ul_us_inx(L, l, ul)
c =unique(L);
NC = length(c);
train = [];
semitrain = [];
test = [];
for ii = 1 : NC
    inx = find(L == c(ii));
    s = sampling(1 : length(inx), l);
    train = [train;inx(s)];
    inx = setdiff(inx, inx(s));
    s = sampling(1 : length(inx), ul);
    semitrain = [semitrain;inx(s)];
    inx = setdiff(inx, inx(s));
    test = [test;inx];
end
end

function s = sampling(R, n)
% 选择抽样，R为记录集合，n为抽取的样本数
% 算法参考：D. E. Knuth, TAOCP, vol.2, pp142，稍有改动
  
% 编写函数时用的测试数据
if ~nargin
    R = 1 : 8;
    n = 4;
end
  
N = length(R);
t = 0;   % 处理过的记录总数
m = 0;   % 已选得的记录数
  
while 1
    U  = rand;
    if (N-t)*U < n-m
        m = m + 1;
        s(m) = R(t+1);
        % 若已抽取到足够的记录，则算法终止
        if m >= n, break, end
    end
    t = t + 1;
end
end

function [Xtrain, Xtest] = dimension_reduction_sub(Xtrain, Ltrain, XULTrain, Xtest, sel_dim, reduce_type)
switch(reduce_type)
    case {'SSMCFS', 'SMCFS','LEASTFS'}
        [ resel_label ] = semi_selected_main(Xtrain, Ltrain, XULTrain, sel_dim, reduce_type  );
        Xtrain = Xtrain(:, resel_label);
        if ~isempty(Xtest)
            Xtest = Xtest(:, resel_label); 
        end
    case 'PCA'
        XPCA = [Xtrain;Xtest];
        options.ReducedDim = sel_dim;
        [eigvector, eigvalue] = PCA(XPCA, options);
        Xtrain = Xtrain * eigvector;
        Xtest = Xtest * eigvector;
    case 'LDA'
        [eigvector, eigvalue] = LDA(Ltrain, [], Xtrain);
        Xtrain = Xtrain * eigvector;
        Xtest = Xtest * eigvector;
    case {'mrmr', 'SFS', 'MCFS', 'LaplacianScore', 'fcbf', 'relief', 'condred', 'cmi', 'disr'}
        [resel_label] = selected_main(Xtrain,Ltrain, sel_dim, reduce_type);
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);     
end
end

%     [L,~] = lmnn2(XLTrain', LTrain',9);
%     XLTrain = XLTrain';
% XLTrain=L*XLTrain;
% XLTrain = XLTrain';
% XTest = XTest';
% XTest = L*XTest;
% XTest = XTest';
%[XLTrain, XULTrain, XTest] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, XTest,'LDA', 200);

%     [ p  ] = locate_face_reco( XLTrain, LTrain, XULTrain, XTest, LTest);
%     pknn = [pknn;p];
%     options.PCARatio = 0.98;
%     [eigvector, eigvalue] = PCA(X, options);
%     XLTrain = XLTrain * eigvector;
%     XULTrain = XULTrain * eigvector;
%     XTest = XTest * eigvector;

