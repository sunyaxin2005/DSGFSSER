function [  ] = feature_stat_locate( features, feature_type, save_dir_main, neibor_num)
%FEATURE_STAT_LOCATE Summary of this function goes here
%   Detailed explanation goes here
[files_num, cluster_num] = size(features);
all_stats_res = [];
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

% function [ C ] = locate_stats_cluster( all_sentences_features )
% %LOCATE_STATS_CLUSTER Summary of this function goes here
% %   Detailed explanation goes here
% files_num = length(all_sentences_features);
% all_features = [];
% for inx = 1 : files_num
%     features = all_sentences_features{inx};
%     all_features = [all_features; features'];
% end
% [IDX,C] = kmeans(all_features, 5);
% end

% function [ all_stat_resss ] = stat_using_locate( features)
% stat_type = 'std';
% process_type = 'origin';
% [frame_num features_dim] = size(features);
% %stat_res = zeros(features_dim, 4);
% all_stat_resss = [];
% for feature_index = 1 : 2: features_dim-3
%     stat_resss = [];
%     mean_stat_resss = [];
%     mean_diff_stat = [];
%     std_stat = [];
%     std_diff_stat = [];
%     for frame_inx = 3 : 1: frame_num-2  
%         sig_features = features(frame_inx -2 : frame_inx + 2, feature_index);
%         sig_features = sig_features(~isnan(sig_features));
%         stat_res1 = spr_stat_one_dim(sig_features, stat_type, process_type);
%         
%         sig_features = features(frame_inx -2 : frame_inx + 2, feature_index + 1);
%         sig_features = sig_features(~isnan(sig_features));
%         stat_res2 = spr_stat_one_dim(sig_features, stat_type, process_type);
%         
%         sig_features = features(frame_inx -2 : frame_inx + 2, feature_index + 2);
%         sig_features = sig_features(~isnan(sig_features));
%         stat_res3 = spr_stat_one_dim(sig_features, stat_type, process_type);
%         
%         stat_res = [stat_res1(:)'; stat_res2(:)';stat_res3(:)'];
%         stat_res = max(stat_res);
%         mean_stat_resss = [mean_stat_resss, stat_res(1)];
%         mean_diff_stat =  [mean_diff_stat, stat_res(3)];
%         std_stat = [std_stat, stat_res(2)];
%         std_diff_stat = [std_diff_stat, stat_res(4)];
%     end
%     dim = size(mean_stat_resss, 2);   
%     %[B, INX] = sort(std_diff_stat, 'descend');
%     [B, INX] = sort(mean_stat_resss, 'descend');
%     %B = mean_stat_resss(INX);
%     stat_resss = [stat_resss, mean(B(1 : floor(dim/5))), mean(B(floor(dim/5) : floor(dim * 2/5))), mean(B(floor(dim*2/5) : floor(dim * 3/5)))];
%     [B, INX] = sort(mean_diff_stat, 'descend');
%     %B = mean_diff_stat(INX);
%     stat_resss = [stat_resss, mean(B(1 : floor(dim/5))), mean(B(floor(dim/5) : floor(dim * 2/5))), mean(B(floor(dim*2/5) : floor(dim * 3/5)))];
%      [B, INX] = sort(std_stat, 'descend');
%     %B = std_stat(INX); 
%     stat_resss = [stat_resss, mean(B(1 : floor(dim/5))), mean(B(floor(dim/5) : floor(dim * 2/5))), mean(B(floor(dim*2/5) : floor(dim * 3/5)))];
%     [B, INX] = sort(std_diff_stat, 'descend');
%     %B = std_diff_stat(INX); 
%     stat_resss = [stat_resss, mean(B(1 : floor(dim/5))), mean(B(floor(dim/5) : floor(dim * 2/5))), mean(B(floor(dim*2/5) : floor(dim * 3/5)))];
%     all_stat_resss = [all_stat_resss;stat_resss];
% end
% 
% end

% function [ all_stat_resss ] = stat_using_locate1( mean_stats,  mean_diff_stats, std_stats, std_diff_stats, mean_stats_C,  mean_diff_stats_C, std_stats_C, std_diff_stats_C)
% [features_dim frame_num ] = size(mean_stats);
% all_stat_resss = [];
% 
% dist_mean = dist2(mean_stats', mean_stats_C);
% [C,I1] = min(dist_mean');
% dist_mean_diff = dist2(mean_diff_stats', mean_diff_stats_C);
% [C,I2] = min(dist_mean_diff');
% dist_std = dist2(std_stats', std_stats_C);
% [C,I3] = min(dist_std');
% dist_std_diff = dist2(std_diff_stats', std_diff_stats_C);
% [C,I4] = min(dist_std_diff');
% 
% for feature_index = 1 : 1: features_dim
%     stat_resss = [];
%     mean_stat = (mean_stats(feature_index , :));
%     mean_diff_stat = (mean_diff_stats(feature_index, :));
%     std_stat = (std_stats(feature_index, :));
%     std_diff_stat = (std_diff_stats(feature_index, :));
%     
%     INX = I1;
%     B = mean_stat;
%     stat_resss = [stat_resss, mean(B(INX == 1)), mean(B(INX == 2)), mean(B(INX == 3)), mean(B(INX == 4)), mean(B(INX == 5))];
%     INX = I2;
%     B = mean_diff_stat;
%     stat_resss = [stat_resss, mean(B(INX == 1)), mean(B(INX == 2)), mean(B(INX == 3)), mean(B(INX == 4)), mean(B(INX == 5))];
%     INX = I3;
%     B = std_stat;
%     stat_resss = [stat_resss, mean(B(INX == 1)), mean(B(INX == 2)), mean(B(INX == 3)), mean(B(INX == 4)), mean(B(INX == 5))];
%     INX = I4;
%     B = std_diff_stat;
%     stat_resss = [stat_resss, mean(B(INX == 1)), mean(B(INX == 2)), mean(B(INX == 3)), mean(B(INX == 4)), mean(B(INX == 5))];
%     
% %     dim = size(mean_stat, 2);   
% %     %[B, INX] = sort(std_diff_stat, 'descend');
% %     
% %     [B, INX] = sort(mean_stat, 'descend');
% %     %B = mean_stat(INX1);
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(mean_diff_stat, 'descend');
% %     %B = mean_diff_stat(INX2);
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(std_stat, 'descend');
% %     %B = std_stat(INX3); 
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(std_diff_stat, 'descend');
% %     %B = std_diff_stat(INX4); 
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
%     all_stat_resss = [all_stat_resss;stat_resss];
% end
% end
