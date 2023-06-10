function [ FeaIndex,FeaNumCandi, cLbest ] = feature_select_tradition_main( A, emotion_label, FeaNumCandi, feature_sel_type, Para )
%FEATURE_SELECT_TRADITION_MAIN Summary of this function goes here
%   Detailed explanation goes here
switch feature_sel_type
    case 'FS'
        if ~isempty(emotion_label)
            [out] = fsFisher(A,emotion_label);
            eigvec =  out.fList;
        else
             W = constructW(A);
             Y = LaplacianScore(A, W);
             eigvec = -Y;
        end
    case 'LLFS'
        if exist('Para', 'var')
            tmpPara.sigma = 1;
            tmpPara.lambda = Para.lambda;
            tmpPara.plotfigure = 0;
        else
            tmpPara.sigma = 1;
            tmpPara.lambda = 0.1;
            tmpPara.plotfigure = 0; 
        end   
        eigvec = Logo_MulticlassProblem(A', emotion_label, tmpPara);
        
    case 'LS'
        W = constructW(A);
        Y = LaplacianScore(A, W);
        eigvec = -Y;
end
[dump,idx] = sort(eigvec,'descend');  
FeaIndex = cell(1,length(FeaNumCandi));
for i = 1:length(FeaNumCandi)      
     FeaIndex{i} = idx(1:FeaNumCandi(i));
end

SelectFeaIdx = FeaIndex{1};  
cLbest = SelectFeaIdx;
end

