function [ stats_res, weights_res] = gmm_supervector( M_speakers)
%GMM_CLASS Summary of this function goes here
%   Detailed explanation goes here
file_num = length(M_speakers);
stats_res = [];
weights_res = [];
for file_inx = 1 : file_num
    %[B,IX] = sort(M_speakers(file_inx).data.PComponents, 'descend');
    stat_res = M_speakers(file_inx).data.mu;
    stat_res = stat_res';
    stats_res = [stats_res; stat_res(:)'];
    weight_res = repmat(M_speakers(file_inx).data.PComponents, 36, 1);
    weights_res = [weights_res;weight_res(:)'];
end
%feature_select( stats_res,  files_emotion_label, save_sel_res );
end

