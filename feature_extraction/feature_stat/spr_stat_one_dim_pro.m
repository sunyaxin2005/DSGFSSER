function [ stat_res, dim_stat_res ] = spr_stat_one_dim_pro(features, sel_stat_dim, speech_label)
%SUN_STAT Summary of this function goes here
%   Detailed explanation goes here
%features, 暂时定义为frame_num * 1, 
%stat_type, 目前有'mean', 'std', 'centroid', 'entropy', 'flatness', 'hue2'
%均值，方差，重心，光滑度，二阶hue矩
%process_type, 目前有'origin','deriv1', 'deriv2', 'deriv3', 'deriv4'
%‘不经过任何预处理’，‘一阶导’，‘二阶导’， ‘三阶导’， ‘四阶导’
if(nargin == 1)
   sel_stat_dim = [1:24];
   %sel_stat_dim = sel_stat_dim + 1;
elseif sel_stat_dim == 1
    sel_stat_dim = [1:24];
end
if(nargin == 2 || nargin == 1)
   useing_label = 0;
else
   useing_label = 1;
end
features_dim = size(features, 2);
stat_res = zeros(features_dim, 24);
dim_stat_res = zeros(features_dim, 24);
for feature_index = 1 : features_dim
    feature = features(:, feature_index);
    %使用label
    if useing_label
       undiff_stats_featue = feature(speech_label == 1); 
    else
       undiff_stats_featue = feature;
    end
    undiff_stats_featue = undiff_stats_featue(~isnan(undiff_stats_featue));
    undiff_stats_featue = undiff_stats_featue(~isinf(undiff_stats_featue));
    if length(undiff_stats_featue) < 3
        stat_res(feature_index, :) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        continue;
    end
    [mean_res std_res max_res min_res k s range_res mid_value] = max_min_mean_std_feature( undiff_stats_featue);
    stat_res(feature_index, 1 : 8) = [mean_res, std_res, max_res, min_res, k, s, range_res, mid_value];
    diff_stats_featue = deriv1( feature);   
    diff_diff_stats_featue = deriv1(diff_stats_featue);
    if useing_label
       diff_stats_featue = diff_stats_featue(speech_label == 1); 
    end
    diff_stats_featue = diff_stats_featue(~isnan(diff_stats_featue));
    diff_stats_featue = diff_stats_featue(~isinf(diff_stats_featue));
    [mean_res std_res max_res min_res k s range_res mid_value] = max_min_mean_std_feature(diff_stats_featue);
    stat_res(feature_index, 9 : 16) = [mean_res, std_res, max_res, min_res, k, s, range_res, mid_value];
     if useing_label
       diff_diff_stats_featue = diff_diff_stats_featue(speech_label == 1); 
    end
    diff_diff_stats_featue = diff_diff_stats_featue(~isnan(diff_diff_stats_featue));
    diff_diff_stats_featue = diff_diff_stats_featue(~isinf(diff_diff_stats_featue));
    [mean_res std_res max_res min_res k s range_res mid_value] = max_min_mean_std_feature(diff_diff_stats_featue);
    stat_res(feature_index, 17 : 24) = [mean_res, std_res, max_res, min_res, k, s, range_res, mid_value];  
    dim_stat_res(feature_index, :) = 1 : 24; 
end
stat_res = stat_res(:, sel_stat_dim);
dim_stat_res = dim_stat_res(:, sel_stat_dim);
end

function [mean_res std_res max_res min_res k s range_res mid_value] = max_min_mean_std_feature( features)
if(isempty(features))
    mean_res = NaN;
    std_res = NaN;
    max_res = NaN;
    min_res = NaN;
    k = NaN;
    s = NaN;
    range_res = NaN;
    mid_value = NaN;
else
    mean_res = mean(features);
    std_res = std(features);   
    max_res = max(features);
    min_res = min(features);
    range_res = max_res - min_res;
    mid_value = median(features);
    [k, s] = this_kurtosis_skewness(features);
end
end
function [ features ] = deriv1( features)
[features_num dim_num] = size(features);
for features_index = 1 : features_num - 1;
    features(features_index, 1) = features(features_index + 1, 1) - features(features_index, 1);
end
features(features_num, 1) = features(features_num - 1, 1);
end
function [k, s] = this_kurtosis_skewness(X)
mean_x = mean(X);
X_U = X - mean_x;
std_x = std(X);
k = mean(X_U.^4)/std_x.^4;
s = mean(X_U.^3)/std_x.^3;
end
% function [max_res min_res mean_res std_res entropy_res] = max_min_mean_feature_diff( features)
%     max_res = 0;
%     min_res = 0;
%     mean_res = 0;
%     std_res = 0;
%     entropy_res = 0;
% if(isempty(features))
%     max_res = 0;
%     min_res = 0;
%     mean_res = 0;
%     std_res = 0;
% else
%     features = deriv1(features);
% %     max_res = max(features);
% %     min_res = min(features);
%     mean_res = mean(features);
%     std_res = std(features);
% %     kurtosis_res = Kurtosis(features);
% %    Skewness_res = Skewness(features);    
% %     entropy_res = spr_entropy(features(:));
% end
% end

% function [ features ] = hue( features)
% [features_num dim_num] = size(features);
% 
% m00 = 0; m01 = 0; m10 = 0; 
% u00 = 0; u02 = 0; u20 = 0;
% n02 = 0; n20 = 0;
% u21 = 0; u12 = 0;
% u30 = 0; u03 = 0;      
% u11 = 0;
%         
% for features_index = 1 : features_num;
%     m00 = m00 + features(features_index, 1);
%     m01 = m01 + features(features_index, 1) * features_index;
% end
% mean_y = m01/m00;
% for features_index = 1 : features_num;
%     u00 = u00 +  features(features_index, 1);
%     u02 = u02 +  features(features_index, 1) * ((row1 - mean_y) * (row1 - mean_y));
% end
% n02 = u02/(u00 * u00);
% end
% function [ Kurtosis_res ] = Kurtosis( features)
% Kurtosis_res = (mean((features - mean(features)))).^4/(std(features)).^4;
% end
% function [ Kurtosis_res ] = Skewness( features)
% Kurtosis_res = (mean((features - mean(features)))).^3/(std(features)).^3;
% end


% function [stat_res] = features_hist(features)
% min_f = min(features);
% max_f = max(features);
% stat_res = zeros(1, 32);
% for ci = 1 : length(features)
%     bin = floor((features(ci) - min_f) * 31/(max_f - min_f)) + 1;
%     stat_res(bin) =  stat_res(bin) + 1;
% end
% [ stat_res ] = spr_entropy(features);
% end

% function [ stat_res ] = spr_mean( features)
% stat_res = mean(features);
% end
% 
% function [ stat_res ] = spr_std( features)
% stat_res = std(features);
% end
% 
% function [ centroid_res ] = spr_centroid( features)
% 
% [features_num cols] = size(features);
% xxx = 0;
% yyy = 0;
% centroid_res = 0;
% for feature_index = 1 : features_num
%    xxx = xxx +  feature_index * features(feature_index, 1)^2;
%    yyy = yyy + features(feature_index, 1)^2;
% end
% if( yyy ~= 0)
%     centroid_res = xxx/yyy;
% end 
% 
% end

% function [ entropy_res ] = spr_entropy( features)
% features = mat2gray(features);
% features = im2double(features);
% entropy_res = 0;
% featuresum = sum(features);
% if featuresum == 0;
%     entropy_res = 0;
%     return;
% end
% pdf_features = features/featuresum;
% [features_num cols] = size(pdf_features);
% for feature_index = 1 : features_num
%     if features(feature_index, 1) <= 0
%         continue;
%     end
%     entropy_res = entropy_res + features(feature_index, 1) * log2(features(feature_index, 1));
% end
% entropy_res = -entropy_res;
% end

