function [ p pc classes] = classifier_main( XLTrain, XTest, LTrain, LTest, classifier_type, XULTrain, ULTrain)
%CLASSIFIER_MAIN Summary of this function goes here
%   Detailed explanation goes here
if nargin < 6
    XULTrain = [];
    ULTrain = [];
end
switch classifier_type
    case {'SVM', 'NN', 'knn', 'fknn', 'hknn', 'LMC', 'lpc', 'alh', 'RLMC', 'CGMC', 'CGMC2'}
        [ p, classes ] = classifier_tradition_main( XLTrain, XTest, LTrain, LTest, classifier_type);
    case {'tradition_src', 'kcckneighbor_src', 'WSRC', 'tradition_src_un_norm',  'WSRC_un_norm', 'CRC_RLS', 'LSDL_SRC', 'WSRC1', 'GCMC_SRC', 'GCMC_SRC0'}
        [ classes] = SRC_MAIN( XLTrain, LTrain, XTest, XULTrain, ULTrain, classifier_type);
        p = length(LTest(classes == LTest))/length(classes);     
    case 'DBN_classifier'
        [ p, bad ] = DBN_classifier( XLTrain, XTest, LTrain, LTest );
        p = 1 - p;
        classes = LTest;
        classes(bad) = -1;
    case 'softmax'
        lambda = 1e-4;
        options.maxIter = 400;
        numClasses = numel(unique(LTrain));
        InputSize=size(XLTrain, 2);
        softmaxModel = softmaxTrain(InputSize, numClasses, lambda, XLTrain', LTrain', options);
        [classes] = softmaxPredict(softmaxModel, XTest');
        classes = classes(:);
        p = length(LTest(classes == LTest))/length(classes);   
    case 'logistic_muticlass'
        [ classes] = logistic_muticlass(XLTrain, LTrain, XTest);
         p = length(LTest(classes == LTest))/length(classes);  
end
UL = unique(LTrain);
pc = zeros(1, length(UL));
for ii = 1 : length(UL)
    classesii = classes(LTest == UL(ii));
    LTestii = LTest(LTest == UL(ii));
    if (isempty(LTestii) ) 
         pc(ii) = 0;
    else
         pc(ii) = length(LTest(LTestii == classesii))/length(classesii);
    end
end
end
% 
% function [ classes] = logistic_muticlass(XLTrain, LTrain, XTest)
% %ÑµÁ·Á½Á½·ÖÀàÆ÷
% UL = unique(LTrain);
% lambda = 1e-4;
% options.maxIter = 400;
% numClasses = 2;
% InputSize=size(XLTrain, 2);
% softmaxModels = cell(length(UL), length(UL));
% for ii = 1 : length(UL)
%     XLTrainii = XLTrain(LTrain ==UL(ii), : );
%     LTrainii = LTrain(LTrain ==UL(ii));
%     LTrainii = LTrainii - LTrainii + 1;
%     for jj = ii + 1 : length(UL)
%         XLTrainjj = XLTrain(LTrain ==UL(jj), : );
%         LTrainjj =  LTrain(LTrain ==UL(jj));
%         LTrainjj = LTrainjj - LTrainjj + 2;
%         XLTrainiijj = [XLTrainii;XLTrainjj];
%         LTrainiijj = [LTrainii;LTrainjj];
%         softmaxModels{ii, jj} = softmaxTrain(InputSize, numClasses, lambda, XLTrainiijj', LTrainiijj', options);
%     end
% end
% 
% %²âÊÔ½×¶Î
% classes = zeros(size(XTest, 1), (length(UL) - 1) * (length(UL) - 2)/2);
% num = 0;
% for ii = 1 : length(UL)
%     for jj = ii + 1 : length(UL)
%         num = num + 1;
%         softmaxModel = softmaxModels{ii, jj} ;
%         [class] = softmaxPredict(softmaxModel, XTest');
%         class1 = class;
%         class1(class == 1) = ii;
%         class1(class == 2) = jj;
%         classes(:, num) = class1;
%     end
% end
% [ classes, max_votes_num ] = vote_by_result( classes' );
% classes = classes';
% end
% 
% function [ vote_res, max_votes_num ] = vote_by_result( reco_result )
% %VOTE_BY_RESULT Summary of this function goes here
% %   Detailed explanation goes here
% [ticket_num, sample_num] = size(reco_result);
% %vote_res = [];
% for samle_inx = 1 : sample_num
%     fea = reco_result(:, samle_inx);
%     u = unique(fea);
%     max_num = 0;
%     max_label = 0;
%     for inx = 1 : length(u)
%         this_num = length (fea(fea == u(inx)));
%         if this_num > max_num
%             max_num = this_num;
%             max_label = u(inx);
%         end
%     end
%     vote_res(samle_inx) = max_label;
%     max_votes_num(samle_inx) = max_num;
% end
% 
% end