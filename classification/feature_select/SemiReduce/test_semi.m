function [ prs, pstds] = test_semi( X, L, ALg_type, sel_dims)
%TEST_SEMI Summary of this function goes here
%   Detailed explanation goes here
%ALg_type = 'LMNN_SWST'
global filename1;
global file_name
ALg_types = {'L21FLDA', 'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'MMP', 'locateLRGPSSDR', 'LRGAAFF', 'LRGPSPSSDA', 'LRGPSSDRCOS'};
random_times = 20;
%X = NormalizeFea(X);

% mind = min(X,[],  2);
% maxd = max(X, [], 2);
% mind = repmat(mind, 1, size(X, 2));
% maxd = repmat(maxd, 1, size(X, 2));
% X = (X - mind)./(maxd - mind);
% X(isnan(X)) = 0;
% X(isinf(X)) = 0;

X = NormalizeFea(X);
stdd = std(X);
meandd = mean(X);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;
X(isnan(X)) = 0;
X(isinf(X)) = 0;
% options.PCARatio = 1.0;
% [eigvector, eigvalue] = PCA(X, options);
% X = X * eigvector;

p1 = zeros(random_times, length(sel_dims) + 1);
p2 = zeros(random_times, length(sel_dims) + 1);
p3 = zeros(random_times, length(sel_dims) + 1);
p4 = zeros(random_times, length(sel_dims) + 1);
p5 = zeros(random_times, length(sel_dims) + 1);
p6 = zeros(random_times, length(sel_dims) + 1);
p7 = zeros(random_times, length(sel_dims) + 1);
p8 = zeros(random_times, length(sel_dims) + 1);
p9 = zeros(random_times, 1);
p11 = zeros(random_times, 1);
p10 = zeros(random_times, length(sel_dims), 100);
for ii = 1 : random_times
    %try
    l = 6;
    ul = 3;
    us = 30;
    [train, semitrain, test] = get_l_ul_us_inx(L, l, ul, us);
    %[train, semitrain, test] = get_l_ul_us_inx2(L, 0.4, 0.1);
    XLTrain = X(train, :);
    LTrain = L(train);
    XULTrain = X(semitrain, :);
    XTest = X(test, :);
    LTest = L(test);
    ULTrain = L(semitrain);
    %[eigvector, W ] = OPSR( XLTrain, LTrain );
   % tic
 [eigvector, W] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, XTest, ALg_type);
    %toc
     [ prediction1,  accuracy1] = LCKSVD( XLTrain,  LTrain, XTest,  LTest);
     p9(ii) = accuracy1;
%[ XLTrain, XTest] = slda_use(XLTrain, LTrain, XTest);
     UL = unique(LTrain);
        K1 = 10000;
        for z = 1 : length(UL)
            if ( K1 > length(LTrain(LTrain == UL(z))))
                K1 = length(LTrain(LTrain == UL(z)));
            end
        end
       %p10 = p10(:, :, 1 : length(K1: 3 * K1)); 
%      [ test_results ] = lrrSRC( XLTrain, LTrain,  XTest, LTest, sel_dims);
%      sum(test_results)/length(LTest);
try
     [ ID ] = FDDL_RECO( XLTrain,  LTrain, XTest);
             p = length(LTest(ID' == LTest))/length(LTest);
        p11(ii) = p;
catch
end
    for jj = 1 : length(sel_dims)

        str = sprintf('%d_%d', ii, jj);
        filename1 = strcat(file_name, str);
        
         sel_dim = sel_dims(jj);
        if sel_dim > size(W, 2)
            sel_dim = size(W, 2);
        end
        XLTrain1 = XLTrain * eigvector * W(:, 1 : sel_dim);
        XULTrain1 = XULTrain * eigvector * W(:, 1 : sel_dim);
        XTest1 = XTest * eigvector * W(:, 1 : sel_dim); 
                
%          XLTrain1 = XLTrain ;
%         XULTrain1 = XULTrain;
%         XTest1 = XTest ;      
%         [XTest1] = GCMDistance(XTest1, XLTrain1, LTrain, K1);
%         [XLTrain1] = GCMDistance(XLTrain1, XLTrain1, LTrain, K1);
%         XTestC =  pdist2(XTest1, [XLTrain1]);
%         XTrainC =  pdist2(XLTrain1,[XLTrain1]);
%         XULTrainC = pdist2([XLTrain1], [XLTrain1]);        
% %         
%         sorted1 = sort(XTestC, 2);
%         sorted1 = sorted1(:, K1*2+1);
%         sorted2 = sort(XTrainC, 2);
%         sorted2 = sorted2(:, K1*2+1);    
%         sorted3 = sort(XULTrainC, 2);
%         sorted3 = sorted3(:, K1*2+1); 
%         
%         XTest1 = exp(-(XTestC.^2)./(sorted1 * sorted3'));
%         XLTrain1 =  exp(-(XTrainC.^2)./(sorted2 * sorted3'));
%         XULTrain1 = exp(-(XULTrainC.^2)./(sorted3 * sorted3'));
       
        [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
        p = length(LTest(classes == LTest))/length(classes);
        p1(ii, jj) = p;
        
%        % K1 = floor(K1/2);
% % 

%         
%                 [disC] = get_nearst_kclass(XLTrain1, LTrain);
%         [disC1] = get_nearst_kclass(XTest1, LTest);
        
%         Para = [];
%         [FeaIndex,FeaNumCandi, cLbest ] = feature_sel_batch_main( XLTrain1, LTrain, size(XLTrain1, 2)/2, 'MCFS', Para );
%         XLTrain1 = XLTrain1(:, FeaIndex{1});
%         XTest1 = XTest1(:, FeaIndex{1});   
%         XULTrain1 = XULTrain1(:, FeaIndex{1});  
%         
       %[ classes ] = LMC1( XTest1', XLTrain1', LTrain, 5);
%         [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
%         p = length(LTest(classes == LTest))/length(classes);
%         p4(ii, jj) = p;
%         [ I1 ] = class2bmp( XLTrain, LTrain );
%         I1 = mat2gray(I1);
%         %figure, imshow(I1);
%         figure, plot(XLTrain(:, 1), XLTrain(:, 2), '*');
%         
%         [ I1 ] = class2bmp( XLTrain1, LTrain );
%         I1 = mat2gray(I1);
%         figure, plot(XLTrain1(:, 1), XLTrain1(:, 2), '*');
        %figure, imshow(I1);
%         model = svmtrain(LTrain,XLTrain1,  '-t linear');
%         [classes, accuracy, dec_values] = svmpredict(LTest, XTest1, model);  
%         p = length(LTest(classes == LTest))/length(classes);
%         p2(ii, jj) = p; 
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'WSRC1');
%         p = length(LTest(classes == LTest))/length(classes);
%         p7(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'kcckneighbor_srctest');
% %         p = length(LTest(classes == LTest))/length(classes);
% %         p8(ii, jj) = p;        
%          p = [];
%         for kk = 1 : size(classes, 2)
%             temp = length(LTest(classes(:, kk) == LTest))/length(classes);
%             p = [p;temp];
%         end
%         p10(ii, jj, :) = p';
%         [ vote_res, max_votes_num ] = vote_by_result( classes' );
%         p1 = [p;length(LTest(vote_res' == LTest))/length(classes)];
%         p8(ii, jj) = max(p);
% %        
        [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'tradition_src');
        p = length(LTest(classes == LTest))/length(classes);
        p3(ii, jj) = p;
% % % % % % 
        [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'kcckneighbor_src');
        p = length(LTest(classes == LTest))/length(classes);
        p8(ii, jj) = p;
% % % %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'CLSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p4(ii, jj) = p;
% % % % % %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'WSRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p5(ii, jj) = p;
% %         
% % % % % %         
% % % % % %         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, 'tradition_src1');
% % % % % %         p = length(LTest(classes == LTest))/length(classes);
% % % % % %         p6(ii, jj) = p;        
% % % % % %         
        [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'LSDL_SRC');
        p = length(LTest(classes == LTest))/length(classes);
        p4(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, ULTrain, 'two_step_sparse');
%         p = length(LTest(classes == LTest))/length(classes);
%         p6(ii, jj) = p;
%         

%         

% %         

%         
        
        if sel_dims(jj) > size(W, 2)
            break
        end
    end
%     catch
%         %ii = ii - 1;
%     end
end
p1s = std(p1, 0, 1);
p2s = std(p2, 0, 1);
p3s = std(p3, 0, 1);
p4s = std(p4, 0, 1);
p5s = std(p5, 0, 1);
p6s = std(p6, 0, 1);
p7s = std(p7, 0, 1);
p8s = std(p8, 0, 1);

p1 = mean(p1, 1);
p2 = mean(p2, 1);
p3 = mean(p3, 1);
p4 = mean(p4, 1);
p5 = mean(p5, 1);
p6 = mean(p6, 1);
p7 = mean(p7, 1);
p8 = mean(p8, 1);
p1(length(sel_dims) + 1) = max(p1);
[~, inxxx] = max(p1);
p1s(length(sel_dims) + 1) = p1s(inxxx);

p2(length(sel_dims) + 1) = max(p2);
[~, inxxx] = max(p2);
p2s(length(sel_dims) + 1) = p2s(inxxx);

p3(length(sel_dims) + 1) = max(p3);
[~, inxxx] = max(p1);
p3s(length(sel_dims) + 1) = p3s(inxxx);

p4(length(sel_dims) + 1) = max(p4);
[~, inxxx] = max(p3);
p4s(length(sel_dims) + 1) = p4s(inxxx);

p6(length(sel_dims) + 1) = max(p6);
[~, inxxx] = max(p1);
p6s(length(sel_dims) + 1) = p6s(inxxx);

p5(length(sel_dims) + 1) = max(p5);
[~, inxxx] = max(p3);
p5s(length(sel_dims) + 1) = p5s(inxxx);

p7(length(sel_dims) + 1) = max(p7);
[~, inxxx] = max(p1);
p7s(length(sel_dims) + 1) = p7s(inxxx);

p8(length(sel_dims) + 1) = max(p8);
[~, inxxx] = max(p8);
p8s(length(sel_dims) + 1) = p8s(inxxx);

prs = [p1; p2; p3; p4; p5; p6;p7;p8];
pstds = [p1s; p2s; p3s; p4s; p5s; p6s;p7s;p8s];
[mean(p9), std(p9)]
[mean(p11), std(p11)]
% p10 = mean(p10, 1);
% p11 = zeros(size(p10, 2), size(p10, 3));
% p11(:, :) = p10(1, :, :)
end

function [train, semitrain, test] = get_l_ul_us_inx(L, l, ul, us)
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
    if length(inx) > us  
        s = sampling(1 : length(inx), ul);
        test = [test;inx(s)];
    else
        test = [test;inx];
    end
end
end
function [train, semitrain, test] = get_l_ul_us_inx2(L, lratio, ulratio)
c =unique(L);
NC = length(c);
train = [];
semitrain = [];
test = [];
for ii = 1 : NC
    inx = find(L == c(ii));
    l = floor(length(inx) * lratio);
    ul = floor(length(inx) * ulratio);
    s = sampling(1 : length(inx), l);
    train = [train;inx(s)];
    inx = setdiff(inx, inx(s));
    s = sampling(1 : length(inx), ul);
    semitrain = [semitrain;inx(s)];
    inx = setdiff(inx, inx(s));
    test = [test;inx];
end
end

function [train, semitrain, test] = get_l_ul_us_inx3(L, l, ul)
c =unique(L);
NC = length(c);
l = floor(l * NC);
ul = floor(ul * NC);
train = [];
semitrain = [];
test = [];
s = randperm(NC);
inx = 1 : length(L);
for ii = 1 : l
    train = [train,inx(L == c(s(ii)))];
end
for ii = l+1 : l+ul
    semitrain = [semitrain,inx(L == c(s(ii)))];
end
for ii = l+ul + 1 : NC
    test = [test,inx(L == c(s(ii)))];
end
end
function [WS]=train_locate_W(XTrain, XULTrain, LTrain, ULTrain, k)
[disC] = get_nearst_kclass(XTrain, LTrain);
[B,IX] = sort(disC, 2);
WS = [];
ul = unique(LTrain);

for ii = 1 : length(ul)
    XTrainii = [];
    LTrainii = [];
    XULTrainii = [];
    for jj  = 1 : k + 1
        XTrainii = [XTrainii;XTrain(LTrain == ul(IX(ii, jj)), :)];
        Lii = zeros(length(LTrain(LTrain == ul(IX(ii, jj)))), 1) + ul(IX(ii, jj));
        LTrainii = [LTrainii;Lii];
        XULTrainii = [XULTrainii;XULTrain(ULTrain == ul(IX(ii, jj)), :)];
    end
    [eigvector, W] = semi_dimension_reduction_main(XTrainii, LTrainii, XULTrainii, 'LRGPSSDR');
    data.eigvector = eigvector;
    data.W = W;
    data.XTrainii = XTrainii;
    data.LTrainii = LTrainii;
    WS = [WS;data];
end    
end

function [disC] = get_nearst_kclass(XTrain, LTrain)
ul = unique(LTrain);
distance = pdist2(XTrain, XTrain);
disC = zeros(length(ul), length(ul));
for ii = 1 : length(ul)
    for jj = 1 : length(ul)
        disC(ii, jj) = sum(sum(distance(LTrain == ul(ii), LTrain == ul(jj))))/(length(LTrain(LTrain == ul(ii))) * length(LTrain(LTrain == ul(jj))));
    end
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

function [CS] = set_insect2train(IX, WS)
CS = [];
for ii = 1 : length(WS)
    Lii = WS(ii).LTrainii;
    uli = unique(Lii);
    IXi = IX(1:length(uli));
    c = intersect(uli, IXi);
    CS = [CS;length(c)];
end
end

% function [Xtrain, Xtest] = dimension_reduction_sub(Xtrain, Ltrain, XULTrain, Xtest, sel_dim, reduce_type)
% switch(reduce_type)
%     case {'SSMCFS', 'SMCFS','LEASTFS'}
%         [ resel_label ] = semi_selected_main(Xtrain, Ltrain, XULTrain, sel_dim, reduce_type  );
%         Xtrain = Xtrain(:, resel_label);
%         if ~isempty(Xtest)
%             Xtest = Xtest(:, resel_label); 
%         end
%     case 'PCA'
%         XPCA = [Xtrain;Xtest];
%         options.ReducedDim = sel_dim;
%         [eigvector, eigvalue] = PCA(XPCA, options);
%         Xtrain = Xtrain * eigvector;
%         Xtest = Xtest * eigvector;
%     case 'LDA'
%         [eigvector, eigvalue] = LDA(Ltrain, [], Xtrain);
%         Xtrain = Xtrain * eigvector;
%         Xtest = Xtest * eigvector;
%     case {'mrmr', 'SFS', 'MCFS', 'LaplacianScore', 'fcbf', 'relief', 'condred', 'cmi', 'disr'}
%         [resel_label] = selected_main(Xtrain,Ltrain, sel_dim, reduce_type);
%         Xtrain = Xtrain(:, resel_label);
%         Xtest = Xtest(:, resel_label);     
% end
% end

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



%         [ D ] = LSDL_train( XLTrain1, LTrain, K1/2 );
%         LD = [];
%         mark_Ps = [];
%         for z = 1 : length(D)
%             mark_Ps = [mark_Ps, D(z).Djj];
%             UZZ = zeros(size(D(z).Djj, 2), 1) + UL(z);
%             LD = [LD;UZZ];
%         end
%         D = mark_Ps';
      
%XXXX = [XLTrain1;XULTrain1];
% beta = 0.8;
% sparsity_func= 'L1';
% epsilon = [];
% num_iters = 100;
% batch_size = size(XXXX, 1);
% Binit = [];
% num_bases =  size(XXXX, 1)/2;
% fname_save = sprintf('../results/sc_%s_b%d_beta%g_%s', sparsity_func, num_bases, beta, datestr(now, 30));	
% [B S stat] = sparse_coding(XXXX', num_bases, beta, sparsity_func, epsilon, num_iters, batch_size, fname_save, Binit);
% D = B';
% Para = [];
% num_bases =  size(XXXX, 1)/2;
% [FeaIndex,FeaNumCandi, cLbest ] = feature_sel_batch_main( XXXX', [], num_bases, 'MCFS', Para );
% D = XXXX(FeaIndex{1}, :);
        
%         XTestC =  pdist2(XTest1, [D]);
%         XTrainC =  pdist2([D], [D]);
%         XULTrainC = pdist2([D], [D]);


%             UL = unique(LTrain);
%         K1 = 10000;
%         for z = 1 : length(UL)
%             if ( K1 > length(LTrain(LTrain == UL(z))))
%                 K1 = length(LTrain(LTrain == UL(z)));
%             end
%          end
%        % K1 = floor(K1/2);
% % 
%         XTestC =  pdist2(XTest, [XLTrain;XULTrain]);
%         XTrainC =  pdist2(XLTrain,[XLTrain;XULTrain]);
%         XULTrainC = pdist2([XULTrain], [XLTrain;XULTrain]);     
%         XULTrainC1 = pdist2([XLTrain;XULTrain], [XLTrain;XULTrain]);   
% %         
%         sorted1 = sort(XTestC, 2);
%         sorted1 = sorted1(:, K1*4+1);
%         sorted2 = sort(XTrainC, 2);
%         sorted2 = sorted2(:, K1*4+1);    
%         sorted3 = sort(XULTrainC, 2);
%         sorted3 = sorted3(:, K1*4+1); 
%         sorted4 = sort(XULTrainC1, 2);
%         sorted4 = sorted4(:, K1*4+1); 
%         
%         XTest = exp(-(XTestC.^2)./(sorted1 * sorted4'));
%         XLTrain =  exp(-(XTrainC.^2)./(sorted2 * sorted4'));
%         XULTrain = exp(-(XULTrainC.^2)./(sorted3 * sorted4'));
        
%                 [disC] = get_nearst_kclass(XLTrain1, LTrain);
%         [disC1] = get_nearst_kclass(XTest1, LTest);


%  classes = [];
%         for kkk = 1 : length(LTest)
%             uuuuuL = unique(LTrain);
%             ind = find(uuuuuL==LTest(kkk));
%             sel_dim1 = 100000000000;
%            if sel_dim1 > size(WS(ind).W, 2)
%              sel_dim1 = size(WS(ind).W, 2);
%            end
%             XTestkkk = XTest(kkk, :) * WS(ind).eigvector * WS(ind).W(:, 1 : sel_dim1);
%             XLTrain2 = WS(ind).XTrainii * WS(ind).eigvector * WS(ind).W(:, 1 : sel_dim1);
%             
%             LTrain2 = WS(ind).LTrainii;
%             [ classeskkk ] = my1nn( XLTrain2, LTrain2, XTestkkk);
%             classes = [classes;classeskkk];
%             if sel_dim1 > size(W, 2)
%                 break
%             end
% 
%         end
%         p = length(LTest(classes == LTest))/length(classes);
%         p1(ii, jj) = p;      

%      [WS]=train_locate_W(XLTrain, XULTrain, LTrain, ULTrain, length(unique(LTrain))/4);
%          classes = [];
%         for kkk = 1 : length(LTest)
%             uuuuuL = unique(LTrain);
%             ind = find(uuuuuL==LTest(kkk));
%             sel_dim1 = sel_dim;
%            if sel_dim1 > size(WS(ind).W, 2)
%              sel_dim1 = size(WS(ind).W, 2);
%            end
%            [ IX ] = myknn4( XLTrain, LTrain, XTest(kkk, :)  );
%            [CS] = set_insect2train(IX, WS);
%            [C,ind1] = max(CS);
%            ind = ind1(1);
%            
%             XTestkkk = XTest(kkk, :) * WS(ind).eigvector * WS(ind).W(:, 1 : sel_dim1);
%             XLTrain2 = WS(ind).XTrainii * WS(ind).eigvector * WS(ind).W(:, 1 : sel_dim1);
%             
%             LTrain2 = WS(ind).LTrainii;
%             
%             [ classeskkk ] = my1nn( XLTrain2, LTrain2, XTestkkk);
% %             if( ind2 ~= ind)
% %                 bnbnbn = 0;
% %             end
%             classes = [classes;classeskkk];
%             if sel_dim1 > size(W, 2)
%                 break
%             end
%         end
%         p = length(LTest(classes == LTest))/length(classes);
%         p2(ii, jj) = p;  