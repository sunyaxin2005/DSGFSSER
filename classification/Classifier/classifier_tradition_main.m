function [ p, classes ] = classifier_tradition_main( XLTrain, XTest, LTrain, LTest, classifier_type)
%CLASSIFIER_MAIN Summary of this function goes here
%   Detailed explanation goes here
K1 = 10000;
UL = unique(LTrain);
for jj = 1 : length(UL)
    if ( K1 > length(LTrain(LTrain == jj)))
        K1 = length(LTrain(LTrain == jj));
    end
end

XLTrain = XLTrain';
XTest = XTest';
switch classifier_type
    case 'SVM'
        model = svmtrain(LTrain, XLTrain',  '-t polynomial');
        [classes, accuracy, dec_values] = svmpredict(LTest, XTest', model);  
        p = length(LTest(classes == LTest))/length(classes);
    case 'NN'
        [ classes ] = my1nn( XLTrain', LTrain, XTest');
        p = length(LTest(classes == LTest))/length(classes);
    case 'knn'
%         [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
%         k=train_params{1};
        k = 30;
        classes = knnclassify(XTest', XLTrain', LTrain, k);
        p = length(LTest(classes == LTest))/length(classes);
    case 'fknn'
        [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        k=train_params{1};
        [classes,memberships] = fknn1(XLTrain', LTrain', XTest',k);
        p = length(LTest(classes == LTest))/length(classes);
    case 'hknn'
        [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        k=train_params{1};
        lambda=train_params{2};
        classes = zeros(length(LTest), 1);
        for t=1:size(XTest,2)
            classes(t) = hknn1(XLTrain,k,LTrain',lambda, XTest(:,t));
        end;
        p = length(LTest(classes == LTest))/length(classes);
    case 'LMC'
        %[ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        k=30;
        classes=LMC1(XTest,XLTrain,LTrain,k);
        p = length(LTest(classes == LTest))/length(classes);   
    case 'lpc'
        [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        k=train_params{1};
        gamm=train_params{2};
        classes = LPC1(XLTrain,LTrain,XTest,k,gamm);  
        p = length(LTest(classes == LTest))/length(classes);   
    case 'alh'
        [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        T=train_params(1);
        k=train_params(2);
        lambda=train_params(3);
        weight = alhweight(Xtrain', Ltrain', T);
        classes = alhclass(Xtest', Xtrain', Ltrain', weight, k, lambda); 
        p = length(LTest(classes == LTest))/length(classes);   
    case 'RLMC'
        %[ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
        %k=train_params{1};
        k=K1;
        classes=RLMC(XTest,XLTrain,LTrain,k);
        p = length(LTest(classes == LTest))/length(classes);     
    case 'CGMC'
%         [ train_params ] = select_paras( XLTrain, LTrain, classifier_type);
%         k=train_params(1);
%         gamm=train_params(2);
%         eps=train_params(3);
        k = K1;
        classes=CGMC1(XTest, XLTrain,LTrain,k);
        p = length(LTest(classes == LTest))/length(classes);
    case 'CGMC2'
        classes=CGMC2(XTest', XLTrain',LTrain);
        p = length(LTest(classes == LTest))/length(classes);
end

end

