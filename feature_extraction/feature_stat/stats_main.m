function [ stats_features , dim_stat_types ] = stats_main( features, stats_type, speech_labels, locate_stats_res)
%STATS_MAIN Summary of this function goes here
%   Detailed explanation goes here
%features:特征，每一行为一个样本
%class_label:样本类别，n*1
%stats_type目前有两种：'common_stat', 'common_stat_label', 'locate_stat', 'locate_stat_label'
%feature_name:特征的名字，保存中间结果使用, 
%save_dir_main:中间结果保存的主要目录
%locates_dir:局部统计的目录,可选，如果提供则直接使用该目录下的locates_dir，否则重新计算，保存在save_dir_main目录下
%返回值：stats_features，每一行为一个样本
    switch stats_type
        case 'common_stat'
            [stats_features, dim_stat_types ] = common_stat(features);
        case 'common_stat_label'
             [stats_features, dim_stat_types] = common_stat_label( features, speech_labels);
        case 'locate_stat'
            [stats_features, dim_stat_types] = locate_stat(features, locate_stats_res);
        case 'locate_stat_label'
            [ stats_features ] = locate_stat(features, locate_stats_res, speech_labels);
    end
end

function [ all_stats_res, dim_stat_types ] = common_stat(features)
all_stats_res = [];
dim_stat_types = [];
for file_inx = 1 : length(features)
    feature = features{file_inx};  
    feature = (feature - mean(mean(feature)))/(std(feature(:)));
    feature(isnan(feature)) = max(max(feature));
    feature(isinf(feature)) = min(min(feature));
     if ~isempty(feature)
        [stat_res, dim_stat_res] = spr_stat_one_dim_pro(feature);
        all_stats_res = [all_stats_res;stat_res(:)'];
        dim_stat_types = [dim_stat_types;dim_stat_res(:)'];
    else
        stat_res = zeros(24, 1);
        all_stats_res = [all_stats_res;stat_res(:)'];
        dim_stat_types = [dim_stat_types;1:24];
    end
end
end

function [ all_stats_res, dim_stat_res ] = common_stat_label(features, speech_labels)
all_stats_res = [];
for file_inx = 1 : length(features)
    feature = features{file_inx}; 
    speech_label = speech_labels{file_inx};
    if ~isempty(feature)
        if(length(speech_label) > size(feature, 1))
            speech_label = speech_label(1:end-2);
            aaaa = length(speech_label) - size(feature, 1);
            speech_label = speech_label(floor(aaaa/2 + 1):end -floor(aaaa/2));
            if length(speech_label) > size(feature, 1)
                speech_label = speech_label(1:size(feature, 1));
            end
        end
        [stat_res, dim_stat_res] =spr_stat_one_dim_pro( feature, 1, speech_label);
        all_stats_res = [all_stats_res;stat_res(:)'];
    else
        stat_res = zeros(4, 1);
        all_stats_res = [all_stats_res;stat_res(:)'];
    end
end
end

function [ all_stats_res, dim_stat_types] = locate_stat(features, locate_stats_res, speech_labels)
%注意这句话要和所有的地方一致，要不然会重新计算
dim_stat_types = [];
if(nargin == 2)
   useing_label = 0;
else
    useing_label =1;
end
all_stats_res = [];

files_num = length(features);
for file_inx = 1 : files_num
    if useing_label
        speech_label =  speech_labels{file_inx};
        feature = locate_stats_res{file_inx};
        %因为局部统计可能失去边沿，导致和speech_label的维度不一样
        if(length(speech_label) > size(feature{1}, 1))
            aaaa = length(speech_label) - size(feature{1}, 1);
            speech_label = speech_label(floor(aaaa/2 + 1):end -floor(aaaa/2));
        end
        [ stat_res ] = stat_using_locate(feature, speech_label); 
    else
         feature = locate_stats_res{file_inx};
        [ stat_res ] = stat_using_locate(feature);
    end
    if ~isempty(stat_res)
        all_stats_res = [all_stats_res;stat_res(:)'];
        %dim_stat_types = [dim_stat_types;1:24];
    else
        stat_res = zeros(size(features{1}, 2) * 36, 1);
        all_stats_res = [all_stats_res;stat_res(:)'];
    end
end
end
function [ all_stat_resss ] = stat_using_locate( locate_stats_res, speech_label)
%使用label与否判断
if(nargin == 1)
   useing_label = 0;
else
    useing_label =1;
end
all_stat_resss = [];
choose = [1:4, 9 : 12, 17 : 20];
for inx = 1 : length(locate_stats_res)
    locate_stat_res = locate_stats_res{inx};
    %locate_stat_res = (locate_stat_res - mean(mean(locate_stat_res)))/(std(locate_stat_res(:)));
    if useing_label
        locate_stat_res = locate_stat_res(:, speech_label == 1);      
    end
    [ stat_resss ] = get_sort_stats(locate_stat_res);
    all_stat_resss = [all_stat_resss, stat_resss];
end
end

% function [ all_stats_res ] = locate_stat(features, neibor_num, locates_dir, speech_labels)
% %注意这句话要和所有的地方一致，要不然会重新计算
% if(nargin == 3)
%    useing_label = 0;
% else
%     useing_label =1;
% end
% try
%     [ mean_stats,  mean_diff_stats, std_stats, std_diff_stats] = load_locate_stats(locates_dir);
% catch
%     [mean_stats,  mean_diff_stats, std_stats, std_diff_stats ] = save_locate_stats(features, locates_dir, neibor_num);
% end
% all_stats_res = [];
% 
% files_num = length(features);
% for file_inx = 1 : files_num
%     if useing_label
%         speech_label =  speech_labels{file_inx};
%         feature = mean_stats{file_inx};
%         %因为局部统计可能失去边沿，导致和speech_label的维度不一样
%         if(length(speech_label) > size(feature, 1))
%             aaaa = length(speech_label) - size(feature, 1);
%             speech_label = speech_label(floor(aaaa/2 + 1):end -floor(aaaa/2));
%         end
%         [ stat_res ] = stat_using_locate( mean_stats{file_inx},  mean_diff_stats{file_inx}, std_stats{file_inx}, std_diff_stats{file_inx}, speech_label); 
%     else
%         [ stat_res ] = stat_using_locate( mean_stats{file_inx},  mean_diff_stats{file_inx}, std_stats{file_inx}, std_diff_stats{file_inx});
%     end
%     if ~isempty(stat_res)
%         all_stats_res = [all_stats_res;stat_res(:)'];
%     else
%         stat_res = zeros(size(features{1}, 2) * 12, 1);
%         all_stats_res = [all_stats_res;stat_res(:)'];
%     end
% end
% end
% 
% function [ all_stat_resss ] = stat_using_locate( mean_stats,  mean_diff_stats, std_stats, std_diff_stats, speech_label)
% mean_stats = mean_stats';
% mean_diff_stats = mean_diff_stats';
% std_stats = std_stats';
% std_diff_stats = std_diff_stats';
% %使用label与否判断
% if(nargin == 4)
%    useing_label = 0;
% else
%     useing_label =1;
% end
% if useing_label
%     mean_stats = mean_stats(:, speech_label == 1);
%     mean_diff_stats = mean_diff_stats(:, speech_label == 1);   
%     std_stats = std_stats(:, speech_label == 1);
%     std_diff_stats = std_diff_stats(:, speech_label == 1);       
% end
% all_stat_resss = [];
% [ stat_resss ] = get_sort_stats(mean_stats);
% all_stat_resss = [all_stat_resss, stat_resss];
% [ stat_resss ] = get_sort_stats(mean_diff_stats);
% all_stat_resss = [all_stat_resss, stat_resss];
% [ stat_resss ] = get_sort_stats(std_stats);
% all_stat_resss = [all_stat_resss, stat_resss];
% [ stat_resss ] = get_sort_stats(std_diff_stats);
% all_stat_resss = [all_stat_resss, stat_resss];
% %目前使用的是对每一个维度都从高到底的排序
% % [features_dim frame_num ] = size(mean_stats);
% % mean_stats = im2double(imresize(mat2gray(mean_stats), [features_dim/2, frame_num]));
% % mean_diff_stats = im2double(imresize(mat2gray(mean_diff_stats), [features_dim/2, frame_num]));
% % std_stats = im2double(imresize(mat2gray(std_stats), [features_dim/2, frame_num]));
% % std_diff_stats = im2double(imresize(mat2gray(std_diff_stats), [features_dim/2, frame_num]));
% %[features_dim frame_num ] = size(mean_stats);
% %n = hist(mean_stats);
% %[IDX,C] =kmeans(n', 5, 'distance', 'correlation');
% % all_stat_resss = [];
% % for feature_index = 1 : 1: features_dim
% %     stat_resss = [];
% %     mean_stat = (mean_stats(feature_index , :));
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% % end
% %     mean_diff_stat = (mean_diff_stats(feature_index, :));
% %     std_stat = (std_stats(feature_index, :));
% %     std_diff_stat = (std_diff_stats(feature_index, :));
% %     
% %     dim = size(mean_stat, 2);   
% %     [B, INX] = sort(mean_stat, 'descend');
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(mean_diff_stat, 'descend');
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(std_stat, 'descend');
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     [B, INX] = sort(std_diff_stat, 'descend');
% %     stat_resss = [stat_resss, mean(B(1 : floor(dim/3))), mean(B(floor(dim/3) : floor(dim * 2/3))), mean(B(floor(dim*2/3) : floor(dim * 3/3)))];
% %     all_stat_resss = [all_stat_resss;stat_resss];
% % end
% end

function [ stat_resss ] = get_sort_stats(X)
[features_dim frame_num ] = size(X);
did_part_num = 3;
stat_resss = [];
for feature_index = 1 : 1: features_dim
    mean_stat = (X(feature_index , :));
    dim = size(mean_stat, 2);  
    [B, INX] = sort(mean_stat, 'descend');
    stat_res = [];
    for inx = 1 : did_part_num
        start_inx = floor((inx - 1) * dim/did_part_num + 1);
        end_inx = floor(inx * dim/did_part_num);
        stat_res = [stat_res, mean(B(start_inx : end_inx))];
    end
    stat_resss = [stat_resss, stat_res];
end
end


% function [ all_stat_resss ] = stat_using_kmean( mean_stats,  mean_diff_stats, std_stats, std_diff_stats, speech_label)
% mean_stats = mean_stats';
% mean_diff_stats = mean_diff_stats';
% std_stats = std_stats';
% std_diff_stats = std_diff_stats';
% %使用label与否判断
% if(nargin == 4)
%    useing_label = 0;
% else
%     useing_label =1;
% end
% if useing_label
%     mean_stats = mean_stats(:, speech_label == 1);
%     mean_diff_stats = mean_diff_stats(:, speech_label == 1);   
%     std_stats = std_stats(:, speech_label == 1);
%     std_diff_stats = std_diff_stats(:, speech_label == 1);       
% end
% 
% %目前使用的是对每一个维度都从高到底的排序
% all_stat_resss = [kmean_stats(mean_stats);kmean_stats(mean_diff_stats);kmean_stats(std_stats);kmean_stats(std_diff_stats)];
% end
% 
% function [statxs] = kmean_stats(X)
% %dim * frame_num
% cluster_num = 6;
% X(isnan(X)) = 0;
% X(isinf(X)) = 0;
% n = hist(X);
% inx = 1;
% while(1)
%     try
%         inx =inx + 1
%         [IDX,C] =kmeans(n', cluster_num, 'distance', 'cosine');
%         break;
%     catch
%         
%     end
% end
% energy_res = sum(X);
% mean_enerys = zeros(cluster_num, 1);
% for inx = 1 :cluster_num
%     mean_enerys(inx) = mean(energy_res(IDX == inx));
% end
% [B,IX] = sort(mean_enerys);
% statxs = [];
% for inx = 1 : cluster_num
%     statx = mean(X(:, IDX == IX(inx)), 2);
%     statxs = [statxs;statx(:)];
% end
% 
% end