function [  FeaIndex,FeaNumCandi, cLbest ] = feature_sel_batch_main( A, emotion_label, FeaNumCandi, sel_type, Para )
%FEATURE_SEL_BATCH_MAIN Summary of this function goes here
%   Detailed explanation goes here
if ~exist('Para', 'var')
    Para.lambda = 0.8;
    Para.beta = 0.5;
end
switch sel_type
    case {'FS', 'LLFS', 'LS'}
        [ FeaIndex,FeaNumCandi, cLbest ] = feature_select_tradition_main( A, emotion_label, FeaNumCandi, sel_type, Para );
    case {'l1ls_featuresignM', 'l1ls_featuresign', 'mexLasso', 'L21RFS', 'LLESPFS', 'l1ls_featuresign_pami', 'SRM_FS', 'SPFS', 'MCFS', 'l1ls_featuresignM1', 'l1ls_featuresignM2'}
        [ FeaIndex,FeaNumCandi, cLbest] = feature_select_sparse_main( A, emotion_label, FeaNumCandi, sel_type, Para );
    case {'condred', 'relief', 'mrmr', 'fcbf', 'disr', 'cmi', 'cmim', 'jmi', 'SFS'}
        [  FeaIndex,FeaNumCandi, cLbest ] = feature_select_info_main( A, emotion_label, FeaNumCandi, sel_type  );
    case {'two_step_fs'}
        [ FeaIndex,FeaNumCandi, cLbest ] = feature_select_two_step_main( A, emotion_label, FeaNumCandi, sel_type);
end
end

