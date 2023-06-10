function [ prs ] = test_feature_selection(  X, L, alg_type, sel_dims )
%TEST_FEATURE_SELECTION Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(X);
random_times = 2;
 %X = NormalizeFea(X);
stdd = std(X);
meandd = mean(X);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;
X(isnan(X)) = 0;
X(isinf(X)) = 0;
nClass = length(unique(L));

pknn = zeros(random_times, length(sel_dims));
ptrasrc = zeros(random_times, length(sel_dims));
pkccsrc = zeros(random_times, length(sel_dims));
pkccsrc1 = zeros(random_times, length(sel_dims));
pkccsrc2 = zeros(random_times, length(sel_dims));
pkccsrc3 = zeros(random_times, length(sel_dims));
for ii = 1 : random_times
    %try
    l = N * 4/5;
    ul = N * 1/20;
    us = 7;
    %[train, semitrain, test] = get_l_ul_us_inx(L, 3, 5);
    [train, semitrain, test] = get_l_ul_us_inx2(L, 0.3, 0.1);
    XLTrain = X(train, :);
    LTrain = L(train);
    XULTrain = X(semitrain, :);
    XTest = X(test, :);
    LTest = L(test);
    options = [];
    options.nUseEigenfunction = nClass; 
    FeaNumCandi =sel_dims;
    %[eigvector, W] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, ALg_type);
    
    for jj = 1 : length(FeaNumCandi)
%         [XLTrain1, XTest1, XULTrain1] = dimension_reduction_sub(XLTrain, LTrain, XULTrain, XTest, FeaNumCandi(jj),alg_type);
%         [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
%          p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc3(ii, jj) = p;
        
       [XLTrain1, XTest1, XULTrain1] = dimension_reduction_sub(XLTrain, LTrain, XULTrain, XTest, FeaNumCandi(jj),'MCFS');
%                  XLTrain1 = Xtrain;
%         XTest1 = Xtest ; 
        [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
         p = length(LTest(classes == LTest))/length(classes);
        pknn(ii, jj) = p;
%         
         [XLTrain2, XTest2] = dimension_reduction_sub(XLTrain1, LTrain, XULTrain1, XTest1, floor(FeaNumCandi(jj)* 6/10), 'mrmr');
%          [ S ] = remove_redendant( XLTrain1, floor(FeaNumCandi(jj)* 8/10) );
         %XLTrain2 = XLTrain1(:, S);
         %XTest2 = XTest1(:, S);
         [ classes ] = my1nn( XLTrain2, LTrain, XTest2);
        p = length(LTest(classes == LTest))/length(classes);
        pkccsrc(ii, jj) = p;
%         
          [XLTrain2, XTest2] = dimension_reduction_sub(XLTrain1, LTrain, XULTrain1, XTest1, floor(FeaNumCandi(jj)* 7/10), 'mrmr');
         [ classes ] = my1nn( XLTrain2, LTrain, XTest2);
        p = length(LTest(classes == LTest))/length(classes);
        ptrasrc(ii, jj) = p;      
        
        [XLTrain2, XTest2] = dimension_reduction_sub(XLTrain1, LTrain, XULTrain1, XTest1, floor(FeaNumCandi(jj)* 8/10), 'mrmr');
        [ classes ] = my1nn( XLTrain2, LTrain, XTest2);
        p = length(LTest(classes == LTest))/length(classes);
        pkccsrc1(ii, jj) = p;    
        
        [XLTrain2, XTest2] = dimension_reduction_sub(XLTrain1, LTrain, XULTrain1, XTest1, floor(FeaNumCandi(jj)* 9/10), 'mrmr');
        [ classes ] = my1nn( XLTrain2, LTrain, XTest2);
        p = length(LTest(classes == LTest))/length(classes);
        pkccsrc2(ii, jj) = p;     
        
        
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_src');
%         p = length(LTest(classes == LTest))/length(classes);
%         ptrasrc(ii, jj) = p;
%         [XLTrain1, XTest1] = dimension_reduction_sub(XLTrain1, LTrain, XULTrain1, XTest1, 90, 'mrmr');
%          XLTrain1 = Xtrain;
%         XTest1 = Xtest ; 
%         [FeaIndex,FeaNumCandi1] = feature_selection_main(X, FeaNumCandi(jj), options, alg_type);
%         SelectFeaIdx = FeaIndex{1};
%         XLTrain1 = XLTrain(:, SelectFeaIdx);
%          XTest1 = XTest(:, SelectFeaIdx) ; 
        
%         [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
%         p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc(ii, jj) = p;
        
%         [ model,  train_param] = svm_train_main( XLTrain1, LTrain,  'poly_nomial' ); 
%         [ classes ] = svm_test_main( XTest1, LTest, XLTrain1, model,   'poly_nomial');
%          p = length(LTest(classes == LTest))/length(classes);
%           pkccsrc1(ii, jj) = p;
%         [ classes ] = SRC_MAIN( XLTrain, LTrain, XTest, 'CLSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc(ii, jj) = p;
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_src');
%         p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc1(ii, jj) = p;
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'kcckneighbor_src');
%         p = length(LTest(classes == LTest))/length(classes);
%         pkccsrc1(ii, jj) = p;
        
%         if sel_dims(jj) > size(W, 2)
%             break; 
%         end
    end
%     catch
%     end
end
pknn = mean(pknn, 1);
pkccsrc = mean(pkccsrc, 1);
ptrasrc = mean(ptrasrc, 1);
pkccsrc1 =  mean(pkccsrc1, 1);
pkccsrc2 =  mean(pkccsrc2, 1);
pkccsrc3 =  mean(pkccsrc3, 1);
prs = [pknn; pkccsrc; ptrasrc; pkccsrc1;pkccsrc2;pkccsrc3];
%prs = [pknn; pkccsrc];
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
    s = sampling(1 : length(inx), ul);
    semitrain = [semitrain;inx(s)];
    inx = setdiff(inx, inx(s));
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
end
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