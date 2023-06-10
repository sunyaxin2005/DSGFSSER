function [ prs, pstds ] = test_feature_selection_sparse(  X, L, alg_type, sel_dims, Para )
%TEST_FEATURE_SELECTION Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(X);
random_times = 5;
%X = NormalizeFea(X);
stdd = std(X);
inx = stdd ~= 0;
stdd = stdd(inx);
meandd = mean(X);
meandd = meandd(inx);
X = X(:, inx);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;
% % max_X = max(X);
% % min_X = min(X);
% % max_X = repmat(max_X, size(X, 1), 1);
% % min_X = repmat(min_X, size(X, 1), 1);
% % X = (X - min_X)./(max_X- min_X);
X(isnan(X)) = 0;
X(isinf(X)) = 0;
% tol=1e-3;
% meanv = mean(X);
% stdv = std(X);
% meanv = repmat(meanv, size(X, 1), 1);
% stdv = repmat(stdv, size(X, 1), 1);
% Xnorm = (X - meanv + tol)./(stdv + tol);
%Xnorm(isnan(X)) = 0;
%X(isinf(X)) = 0;
% maxv = max(X);
% minv = min(X);
% maxv = repmat(maxv, size(X, 1), 1);
% minv = repmat(minv, size(X, 1), 1);
% X  = (X - minv + tol)./(maxv-minv+tol);


X(isnan(X)) = 0;
X(isinf(X)) = 0;
nClass = length(unique(L));


p1 = zeros(random_times, length(sel_dims) + 1);
p2 = zeros(random_times, length(sel_dims) + 1);
p3 = zeros(random_times, length(sel_dims) + 1);
p4 = zeros(random_times, length(sel_dims) + 1);
p5 = zeros(random_times, length(sel_dims) + 1);
p6 = zeros(random_times, length(sel_dims) + 1);
p7 = zeros(random_times, length(sel_dims) + 1);
p8 = zeros(random_times, length(sel_dims) + 1);
p9 = zeros(random_times, length(sel_dims) + 1);
for ii = 1 : random_times
    %try
    l =3;
    ul = 3;
    us = 7;
    %[train, semitrain, test] = get_l_ul_us_inx(L, l, ul);
    [train, semitrain, test] = get_l_ul_us_inx2(L, 0.8, 0);
    %[train, semitrain, test] = get_l_ul_us_inx3(L, 0.5, 0.0);
    XLTrain = X(train, :);
    LTrain = L(train);
    XULTrain = X(semitrain, :);
    XTest = X(test, :);
    LTest = L(test);
    options = [];
    options.nUseEigenfunction = nClass; 
    FeaNumCandi =sel_dims;
    %[ FeaIndex,FeaNumCandi, cLbest] = feature_select_sparse_main( XLTrain, LTrain, FeaNumCandi, alg_type );
    FDRs = FDR(XLTrain, LTrain );
    XLTrain = XLTrain(:, FDRs > 0.10);
    XTest = XTest(:, FDRs > 0.10);
    [FeaIndex,FeaNumCandi, cLbest ] = feature_sel_batch_main( XLTrain, LTrain, FeaNumCandi, alg_type, Para );
    [ FeaIndexttwo_two ] = train_two_two_sel( XLTrain,  LTrain, FeaNumCandi, alg_type, Para );
    for jj = 1 : length(FeaNumCandi)
%         
        XLTrain1 = XLTrain(:, FeaIndex{jj});
        XTest1 = XTest(:, FeaIndex{jj});   
        XULTrain1 = XULTrain(:, FeaIndex{jj});  

%         XLTrain1 = XLTrain;
%         XTest1 = XTest;   
%         XULTrain1 = XULTrain;   
        
     %   for zzz = 1 : 50
%             [train, semitrain, test] = get_l_ul_us_inx(LTest, l, ul);
%             XLTrain1 = XTest1(train, :);
%             LTrain1 = LTest(train, :);
%             XTest2 = XTest1(test, :);
%             LTest1 = LTest(test, :);
%         XLTrain1 = XLTrain;
%         XTest1 = XTest;   
%         XULTrain1 = XULTrain;   
        
%             [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
%             p = length(LTest(classes == LTest))/length(classes);
%             p1(ii, jj) = p1(ii, jj)+ p;
% [XTest1] = SPPDDistance(XTest1, XLTrain1, LTrain);
% [XLTrain1] = SPPDDistance(XLTrain1, XLTrain1, LTrain);
            [classes1]=LMC1(XTest1',XLTrain1',LTrain,6);
            p = length(LTest(classes1 == LTest))/length(LTest);
            p1(ii, jj) = p1(ii, jj)+ p;
%         [ I1 ] = class2bmp( XLTrain, LTrain );
%         I1 = mat2gray(I1);
%         %figure, imshow(I1);
%         figure, plot(XLTrain(:, 1), XLTrain(:, 2), '*');
%         
%         [ I1 ] = class2bmp( XLTrain1, LTrain );
%         I1 = mat2gray(I1);
%         figure, plot(XLTrain1(:, 1), XLTrain1(:, 2), '*');
        %figure, imshow(I1);
        model = svmtrain(LTrain,XLTrain1,  '-t polynomial');
        [classes2, accuracy, dec_values] = svmpredict(LTest, XTest1, model);  
        p = length(LTest(classes2 == LTest))/length(LTest);
        p2(ii, jj) =  p2(ii, jj)+ p; 
        
        p = length(LTest(classes2 == LTest | classes1 == LTest))/length(LTest);
        p3(ii, jj) =  p3(ii, jj)+ p; 
        
        for kk = 1 : length(LTest)
            if classes1(kk) ~= classes2(kk)
               xtestkk = XTest(kk, :);
               XTrainiijj = FeaIndexttwo_two{classes1(kk), classes2(kk)}.XTrainiijj;
                LTrainiijj = FeaIndexttwo_two{classes1(kk), classes2(kk)}.LTrainiijj;
           FeaIndex1 = FeaIndexttwo_two{classes1(kk), classes2(kk)}.FeaIndex ;
           XLTrain1 = XTrainiijj(:, FeaIndex1{jj});
           XTest1 = xtestkk(:, FeaIndex1{jj});   
           
           [classes]=LMC1(XTest1',XLTrain1',LTrainiijj,6);
           classes1(kk) = classes;

          model = svmtrain(LTrainiijj,XLTrain1,  '-t polynomial');
         [classes, accuracy, dec_values] = svmpredict(1, XTest1, model);  
          classes2(kk) = classes;
       end
        end
       p = length(LTest(classes1 == LTest))/length(LTest);
    p4(ii, jj) = p4(ii, jj)+ p;
   
    p = length(LTest(classes2 == LTest))/length(LTest);
    p5(ii, jj) = p5(ii, jj)+ p;
%         end
%          p1(ii, jj) =  p1(ii, jj)/50;
%          p2(ii, jj) =  p2(ii, jj)/50;
%         p1(ii, jj)
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, XULTrain1, 'tradition_src');
%         p = length(LTest(classes == LTest))/length(classes);
%         p3(ii, jj) = p;
        
%         model = svmtrain(LTrain,XLTrain1,  '-t linear');
%         [classes, accuracy, dec_values] = svmpredict(LTest, XTest1, model);  
%         p = length(LTest(classes == LTest))/length(classes);
%         p3(ii, jj) =  p3(ii, jj)+ p; 
% 
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'LSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p3(ii, jj) = p;
% %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'CLSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p4(ii, jj) = p;
% %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'WSRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p5(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'LSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p6(ii, jj) = p;
        
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_srcLASSO');
%         p = length(LTest(classes == LTest))/length(classes);
%         p7(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'kcckneighbor_srcLASSO');
%         p = length(LTest(classes == LTest))/length(classes);
%         p8(ii, jj) = p;
        
        
%         if sel_dims(jj) > size(W, 2)
%             break
%         end
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
p9s = std(p9, 0, 1);

p1 = mean(p1, 1);
p2 = mean(p2, 1);
p3 = mean(p3, 1);
p4 = mean(p4, 1);
p5 = mean(p5, 1);
p6 = mean(p6, 1);
p7 = mean(p7, 1);
p8 = mean(p8, 1);
p9 = mean(p9, 1);

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
p9(length(sel_dims) + 1) = max(p9);
[~, inxxx] = max(p9);
p9s(length(sel_dims) + 1) = p9s(inxxx);

prs = [p1; p2; p3; p4; p5; p6;p7;p8;p9];
pstds = [p1s; p2s; p3s; p4s; p5s; p6s;p7s;p8s;p9s];
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
    if ul ~= 0;
        s = sampling(1 : length(inx), ul);
        semitrain = [semitrain;inx(s)];
        inx = setdiff(inx, inx(s));
    end

    test = [test;inx];
end
end

function [train, semitrain, test] = get_l_ul_us_inx1(L, l, ul)
inx = 1 : length(L);
    s = sampling(1 : length(L), l);
    train = inx(s);
        inx = setdiff(inx, inx(s));
    s = sampling(1 : length(inx), ul);
    semitrain = inx(s);
        inx = setdiff(inx, inx(s));
    test = inx;
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

function [Xtrain, Xtest, XULTrain] = dimension_reduction_sub(Xtrain, Ltrain, XULTrain, Xtest, sel_dim, reduce_type)
switch(reduce_type)
    case {'SSMCFS', 'SMCFS','LEASTFS'}
        [ resel_label ] = semi_selected_main(Xtrain, Ltrain, XULTrain, sel_dim, reduce_type  );
        Xtrain = Xtrain(:, resel_label);
        if ~isempty(Xtest)
            Xtest = Xtest(:, resel_label); 
        end
        XULTrain = XULTrain(:, resel_label);
    case 'PCA'
        XPCA = [Xtrain;Xtest];
        options.ReducedDim = sel_dim;
        [eigvector, eigvalue] = PCA(XPCA, options);
        Xtrain = Xtrain * eigvector;
        Xtest = Xtest * eigvector;
        XULTrain = XULTrain * eigvector;
    case 'LDA'
        [eigvector, eigvalue] = LDA(Ltrain, [], Xtrain);
        Xtrain = Xtrain * eigvector;
        Xtest = Xtest * eigvector;
        XULTrain = XULTrain * eigvector;
    case {'mrmr', 'SFS', 'MCFS', 'LaplacianScore', 'fcbf', 'relief', 'condred', 'cmi', 'disr', 'jmi', 'cmim', 'MCFSM'}
        [resel_label] = selected_main(Xtrain,Ltrain, sel_dim, reduce_type);
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'SLDFS'
         [ resel_label ] = SLDFS( Xtrain, Ltrain, XULTrain, sel_dim);
         Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'SemiFS'
        [ resel_label ] = SemiFS( Xtrain, Ltrain, XULTrain, sel_dim);
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'NaiveLS'
        [ resel_label ] = NaiveLS( Xtrain, 7, sel_dim );
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'SLS'
        [ resel_label ] = SLS(  Xtrain, Ltrain, 3, sel_dim );
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'LSSFS'
        [ resel_label ] = LSSFS( Xtrain, Ltrain, XULTrain, sel_dim );
         Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case {'l1ls_featuresignM',  'l1ls_featuresign', 'mexLasso', 'SPFS', 'LLESPFS', 'SRM_FS', 'L21RFS', 'l1ls_featuresign_pami'}
        [ resel_label ] = feature_select_sparse_main( Xtrain, Ltrain, sel_dim, reduce_type );
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'SPFS_SFS'%On Similarity Preserving Feature Selection, IEEE TRANSACTIONS ON KNOWLEDGE AND DATA ENGINEERING, 2013
        [ resel_label ] = SPFS_SFS( Xtrain, Ltrain, sel_dim );
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'FS'
        [out] = fsFisher(Xtrain',Ltrain);
        resel_label =  out.fList(1 : sel_dim);
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);
    case 'LLFS'
        Para.sigma = 1;
        Para.lambda = 0.5;
        Para.plotfigure = 0;
        Weight = Logo_MulticlassProblem(Xtrain', Ltrain, Para);
        [~, out.fList] = sort(Weight, 'descend');
         resel_label =  out.fList(1 : sel_dim);
        Xtrain = Xtrain(:, resel_label);
        Xtest = Xtest(:, resel_label);   
        XULTrain = XULTrain(:, resel_label);       
end
end
function [ FDRs ] = FDR( X, L )
%FDR Summary of this function goes here
%   Detailed explanation goes here
UL = unique(L);
mean_ls = [];
std_ls = [];
for inx = 1 : length(UL)
    mean_ls = [mean_ls;mean(X(L == UL(inx), :))];
    std_ls = [std_ls;std(X(L == UL(inx), :)).^2];
end
FDRs = 0;
for inx1 = 1 : length(UL)
    for inx2 = inx1:length(UL)
        if (std_ls(inx1, :) + std_ls(inx2, :)) == 0
            FDRs = 0;
        else
            FDRs = FDRs + (mean_ls(inx1, :) - mean_ls(inx2, :)).^2./ (std_ls(inx1, :) + std_ls(inx2, :));
        end
        
    end
end
FDRs = FDRs * 2/((length(UL) - 1) * length(UL));
end


% 
%          for ii = 1 : size(Y, 2)
%              y = Y(:, ii);
%              X = [fea, y];
%             distance = pdist2(X', X');
%             [sorted,index] = sort(distance);
%             sita = sorted(20 * 2+1, :);
%             sita = sita' *  sita;
%             X = exp(-(distance.^2)./sita);
%             fea1 = X(1: size(X, 1) - 1, :)';
%             y = X(size(X, 1), :)';
%             [ M ] = Get_Weight_Near( fea1, 100);
%             eigvec2 = l1ls_featuresignM (fea1, y, M, 0.6, 0.1);
%             eigvec = [eigvec, eigvec2];
%          end