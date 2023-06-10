function [ locate_stat_ress ] = feature_locate_stats( features, neibor_num )
%FEATURE_LOCATE_STATS Summary of this function goes here
%   Detailed explanation goes here
[files_num, cluster_num] = size(features);
locate_stat_ress = cell(files_num, 1);
for files_index = 1 : files_num
    sprintf('parper_locate%d', files_index)
    feature = features{files_index};
    [frame_num features_dim] = size(feature);
    temp_feature = randi(1000, [100,1]);
    stat_res = spr_stat_one_dim_pro(temp_feature);
    stat_length = length(stat_res);
    locate_stats_res = cell(stat_length, 1);
    for feature_index = 1 : 1: features_dim
        locate_stat_res = [];
        for frame_inx = (neibor_num + 1) : 1: frame_num-neibor_num  
            sig_features = feature(frame_inx -neibor_num : frame_inx + neibor_num, feature_index);
            stat_res = spr_stat_one_dim_pro(sig_features);
            locate_stat_res = [locate_stat_res, stat_res(:)];
        end
        for stat_inx = 1 : stat_length
            locate_stats_res{stat_inx} = [locate_stats_res{stat_inx};locate_stat_res(stat_inx, :)];
        end
    end
    locate_stat_ress{files_index} = locate_stats_res;
end
end

