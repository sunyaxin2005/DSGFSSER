function [ features ] = GMM_UBM_PAPER_DATA( features, gmm_type, speech_labels )
%GMM_UBM_PAPER_DATA Summary of this function goes here
%   Detailed explanation goes here
switch gmm_type
   case 'gmm_ubm_svm_locate'
       %[ features ] = feature_locate_stats(features, locates_dir, neibor_num);
   case 'gmm_ubm_svm_label'
       [ features ] = add_labels(features, speech_labels);
   case 'gmm_ubm_svm_locate_label'
        %[ features ] = feature_locate_stats(features, locates_dir, neibor_num);
        [ features ] = add_labels(features, speech_labels);
    case 'gmm_ubm_svm_diff'
        [ features ] = feature_add_diff(features);
    case 'gmm_ubm_svm_diff_label'
        [ features ] = feature_add_diff(features);
        [ features ] = add_labels(features, speech_labels);
    case 'gmm_ubm_svm_diff_diff'
        [ features ] = feature_add_diff_diff(features, speech_labels);
end

end

%为diff准备数据
function [ diff_features ] = feature_add_diff(features)
diff_features = cell(length(features), 1);
for file_num = 1 : length(features)
    file_features = features{file_num};
    features_dim = size(file_features, 2);
    diff_feature = [];
    for features_inx = 1 : features_dim
        feature = file_features(:, features_inx);
        temp = deriv1( feature);
        diff_feature = [diff_feature, temp];
    end
    diff_features{file_num} = [file_features, diff_feature];
end
end

%为diff_diff准备数据
function [ diff_features ] = feature_add_diff_diff(features, speech_labels)
diff_features = cell(length(features), 1);
for file_num = 1 : length(features)
    file_features = features{file_num};
    file_features(isinf(file_features)) = NaN;
    file_features(isnan(file_features)) = max(max(file_features));
    file_features = (file_features - min(min(file_features)))/(max(max(file_features)) - min(min(file_features)));
    features_dim = size(file_features, 2);
    diff_feature = [];
    diff_diff_feature = [];    
    speech_label = speech_labels{file_num};
    if(length(speech_label) > size(file_features, 1))
       speech_label = speech_label(1:end-2);
       aaaa = length(speech_label) - size(file_features, 1);
       speech_label = speech_label(floor(aaaa/2 + 1):end -floor(aaaa/2));
       if length(speech_label) > size(file_features, 1)
          speech_label = speech_label(1:size(file_features, 1));
        end
    end
    for features_inx = 1 : features_dim
        feature = file_features(:, features_inx);
        temp = deriv1( feature);
        temp1 =  deriv1( temp);
        diff_feature = [diff_feature, temp];
        diff_diff_feature = [diff_diff_feature, temp1];
    end
    file_features = file_features(speech_label == 1, :);
    diff_feature = diff_feature(speech_label == 1, :);
    diff_diff_feature = diff_diff_feature(speech_label == 1, :);
    diff_features{file_num} = [file_features, diff_feature, diff_diff_feature];
end
end

%%%%%%%%%以下函数为局部统计准备数据
function [ locate_stats_features ] = feature_locate_stats(features, locates_dir, neibor_num)
try
[ mean_stats,  mean_diff_stats, std_stats, std_diff_stats] = load_locate_stats(locates_dir);
catch
 [mean_stats,  mean_diff_stats, std_stats, std_diff_stats ] = save_locate_stats(features, locates_dir, neibor_num);
end
file_num = length(mean_stats);
locate_stats_features = cell(file_num, 1);
for file_inx = 1 : file_num
    locate_stats_features{file_inx} = [mean_stats{file_inx}, mean_diff_stats{file_inx}, std_stats{file_inx}, std_diff_stats{file_inx}];
end
end

function [ feature ] = deriv1( feature)
sample_dim = length(feature);
for sample_dim_index = 1 : sample_dim - 1;
    feature(sample_dim_index) = feature(sample_dim_index + 1) - feature(sample_dim_index);
end
feature(sample_dim_index) = feature(sample_dim_index - 1);
end