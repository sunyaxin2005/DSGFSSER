function [ all_stat_resss ] = feature_stats_locate(feature_type, save_dir_main, neibor_num)
%FEATURE_STATS_LOCATE Summary of this function goes here
%   Detailed explanation goes here
str = sprintf('%d', neibor_num);
save_feature_path = strcat(save_dir_main, feature_type, '\', 'locate_stats', str, '\');
if ~exist(save_feature_path, 'dir')
    mkdir(save_feature_path);
end
try
[ mean_stats,  mean_diff_stats, std_stats, std_diff_stats] = load_locate_stats(save_feature_path);
catch
 [mean_stats,  mean_diff_stats, std_stats, std_diff_stats ] = save_locate_stats(features, save_feature_path, neibor_num);
end
[ all_stat_resss ] = stat_using_locate( mean_stats,  mean_diff_stats, std_stats, std_diff_stats );
end

function [ all_stat_resss ] = stat_using_locate( mean_stats,  mean_diff_stats, std_stats, std_diff_stats )
[features_dim frame_num ] = size(mean_stats);
all_stat_resss = [];
for feature_index = 1 : 1: features_dim
    stat_resss = [];
    mean_stat = (mean_stats(feature_index , :));
    mean_diff_stat = (mean_diff_stats(feature_index, :));
    std_stat = (std_stats(feature_index, :));
    std_diff_stat = (std_diff_stats(feature_index, :));
    
    dim = size(mean_stat, 2);   
    %[B, INX] = sort(std_diff_stat, 'descend');
    
    [B, INX] = sort(mean_stat, 'descend');
    %B = mean_stat(INX);
    stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
    [B, INX] = sort(mean_diff_stat, 'descend');
    %B = mean_diff_stat(INX);
    stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
    [B, INX] = sort(std_stat, 'descend');
    %B = std_stat(INX); 
    stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
    [B, INX] = sort(std_diff_stat, 'descend');
    %B = std_diff_stat(INX); 
    stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
    all_stat_resss = [all_stat_resss;stat_resss];
end
end
