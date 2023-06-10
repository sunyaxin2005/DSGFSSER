function [ FeaIndex,FeaNumCandi, cLbest ] = feature_select_two_step_main( A, emotion_label, FeaNumCandi, sel_type  )
%FEATURE_SELECT_TWO_STEP_MAIN Summary of this function goes here
%   Detailed explanation goes here
Para.lambda = 0.8;
Para.beta = 0.5;
[ FeaIndex1, FeaNumCandi1, cLbest1] = feature_select_sparse_main( A, emotion_label, FeaNumCandi * 2, 'MCFS', Para );
for ii = 1 : length(FeaIndex1)
    X = A (:, FeaIndex1{ii});
    inx = FeaIndex1{ii};
    [FeaIndex, FeaNumCandi, cLbest ] = feature_select_info_main( X, emotion_label, FeaNumCandi1(ii), 'mrmr');
    FeaIndex1{ii} = inx(FeaIndex{1});
end
FeaIndex = FeaIndex1;
end

