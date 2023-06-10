function [  FeaIndex,FeaNumCandi, cLbest ] = feature_select_info_main( A, emotion_label, FeaNumCandi, sel_type  )
%FEATURE_SELECT_INFO_MAIN Summary of this function goes here
%   Detailed explanation goes here
FeaIndex = cell(1,length(FeaNumCandi));
% for ii = 1: length(FeaNumCandi)
     sel_dim = FeaNumCandi(length(FeaNumCandi));
    switch sel_type
        case 'condred'
            [cLbest] = feast('condred', sel_dim, A, emotion_label);
        case 'relief'
            [cLbest] = feast('relief', sel_dim, A, emotion_label);
        case 'mrmr'
            [cLbest] = feast('mrmr', sel_dim, A, emotion_label);
        case 'fcbf'
            [cLbest] = feast('fcbf', sel_dim, A, emotion_label, 0.0);
        case 'disr'
            [cLbest] = feast('disr', sel_dim, A, emotion_label);    
        case 'cmi'
            [cLbest] = feast('cmi', sel_dim, A, emotion_label);     
        case 'cmim'
            [cLbest] = feast('cmim', sel_dim, A, emotion_label);
        case 'jmi'
            [cLbest] = feast('jmi', sel_dim, A, emotion_label);        
        case 'SFS'
            cLbest = mc_SequentialForwardSelection(A', emotion_label, 'mc_scatter_matrices', sel_dim);   
    end
    for i = 1:length(FeaNumCandi)      
        FeaIndex{i} = cLbest(1:FeaNumCandi(i));
    end

SelectFeaIdx = FeaIndex{1};  
cLbest = SelectFeaIdx;
%     FeaIndex{ii} = cLbest;
% end

end

