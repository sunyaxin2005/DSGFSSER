function [ output_args ] = feature_extract_main( wave_files_names, wave_files_path)
%FEATURE_EXTRACT_BYPARAMES Summary of this function goes here
%���ļ������ж�������������ʶ����
%wave_files_names:ʹ�õ������ļ�����
%wave_files_path�������ļ����ڵ�·��
%algorithm_type��ʹ�õ��㷨���ͣ�Ŀǰ��'common_stat', 'common_stat_label', 'locate_stat', 'locate_stat_label', 'gmm_ubm_svm', 'gmm_ubm_svm_locate'
%save_dir_main = 'F:\��һƪ����\berlin15\';
save_dir_main = 'F:\��8ƪ����ʵ����\EmoDBtest\';% The save path of extracted features
%save_dir_main = 'F:\��һƪ����\chinese16\';
[ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_benlin_label( wave_files_names, wave_files_path ); % if you want to extract features for EmoDB, please use uncomment this sentence 
%[ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_label_datas( wave_files_names, wave_files_path );% if you want to extract features for SAVEE, please use uncomment this sentence 
%[ files_name files_people_label files_emotion_label files_sentence_label words_mark] = casia_load_label( wave_files_names, wave_files_path );% if you want to extract features for CASIA, please use uncomment this sentence 
%�õ�ѵ������label
% save files_sentence_label files_sentence_label
% cross_type = 'LOSO';
%train_inxs = get_train_inxs(cross_type, files_emotion_label, files_people_label); 

%save files_people_label files_people_label
%feature_type = {'best_pitch', 'plp', 'RASTA-PLP', 'fromant', 'best_pitch', 'fromant'};
%feature_type = {'MFCC', 'LFPC', 'LPCC',  'plp', 'RASTA-PLP' };%�ڶ�ƪ����ʹ�õ�����
feature_type = {'sun_mfcc', 'Q', 'TQ', 'LPCC', 'plp', 'RASTA-PLP'};%
issave = [1, 0, 1, 1, 0, 1];%feature group save flag, 'Q', 'TQ'are in the same group, 'plp', 'RASTA-PLP' are in the same group
feature_param = [60, 80, 80, 40, 40, 40,]; %the best Mel filters M for 'MFCC', 'Q', 'TQ', the best params for 'LPCC', 'plp', 'RASTA-PLP' are set at the place that calling them
bw = [50, 60, 60, 40, 40, 40]; % the best frequency resolution of STFT b for 'MFCC', 'Q', 'TQ'
bws = [60];% the frequency resolution of STFT b, when it is used for parameter settings, it should be set the  parameters you want to test
feature_params = [120]; % the number of Mel filters M, when it is used for parameter settings, it should be set the  parameters you want to test
%feature_param = [40 40, 10, 40, 40, 40, 40, 40, 40, 40, 40, 40];
% sel_stats_dim = {[1, 2, 5, 10, 12, 14],...%MFCC
%                   [2, 3, 5, 6, 7, 8, 16, 18, 22, 24], ...%'PLP'
%                   [2, 4, 6,8, 9, 18, 20], ...%'LPCC'
%                   [1, 7, 8, 9, 12, 16, 18, 24],...%'RASTA-PLP'
%                   [2, 5, 6, 12, 16, 22, 24]};%'hu_dct'
% sel_stats_dim = {[1:24],...%MFCC
%                   [1:24], ...%'PLP'
%                   [1:24], ...%'LPCC'
%                   [1:24],...%'RASTA-PLP'
%                   [1:24]};%'hu_dct'
%load speech_labels
speech_labels = [];
%algorithm_types = {'gmm_ubm_svm', 'gmm_ubm_svm_label', 'gmm_ubm_svm_locate', 'gmm_ubm_svm_locate_label', 'gmm_ubm_svm_diff', 'gmm_ubm_svm_diff_label', 'gmm_ubm_svm_diff_diff'};
algorithm_types = {'common_stat'};
    for ii = 1 : length(feature_params)
        for jj = 1 : length(bws)
for algorithm_inx = 1:length(algorithm_types)%�����������ڲ�ͬͳ����ʱ�Ľ��ʹ��
%    try
    neibor_num = 5;
    save_sel_res = strcat(save_dir_main, 'main_reco_res1\', algorithm_types{algorithm_inx}, '\');
    %save_all_features_path = strcat(save_sel_res, 'all_features.mat');
    if ~exist(save_sel_res, 'dir')
        mkdir(save_sel_res);%�����������ڲ�ͬͳ����ʱ�Ľ��ʹ��
    end
    fprintf('\n%s', algorithm_types{algorithm_inx});
    %all_cross_feas = cell(length(train_inxs), 1);
    all_cross_feas = [];
    dim_stat_typess = [];
    stats_dims = [];
    %str1 = sprintf('M%dB%d',  feature_params(ii), bws(jj));
    str1 = sprintf('M%dB%d',  feature_params(ii), bws(jj));
    for inx =1 : length(feature_type)
        algorithm_type = algorithm_types{algorithm_inx};
        fprintf(' %s ', feature_type{inx});
        %����Ŀ¼
        save_dir = strcat(save_dir_main, feature_type{inx}, '\');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        %��ͬ�㷨Ŀ¼
        save_algo_dir = strcat(save_dir, algorithm_types{algorithm_inx}, '\');
        if ~exist(save_algo_dir, 'dir')
            mkdir(save_algo_dir);
        end   
        %��������·��
        save_feature_path = strcat(save_dir, feature_type{inx}, 'features', '.mat');
        %�ֲ�ͳ�Ƶ��м�����ŵ�·��
        str = sprintf('locate_stats%d', neibor_num);
        save_locates_path = strcat(save_dir, str, '.mat');
         %�ֲ������ŵ�·��
        save_stats_path = strcat(save_dir, feature_type{inx}, algorithm_types{algorithm_inx}, str1, '.mat');
        %������ȡ
        try 
            load(save_feature_path);
        catch      
            %feature extraction
            %[ features ] = feature_paper( files_name, feature_type{inx}, feature_params(ii), bws(jj));
            [ features ] = feature_paper( files_name, feature_type{inx}, feature_param(inx), bw(inx));
            save(save_feature_path, 'features');
        end
        locate_stats_res = [];
        if strcmp('locate_stat', algorithm_type)
            try
                load(save_locates_path);
            catch
                [ locate_stats_res ] = feature_locate_stats( features, neibor_num );
                save(save_locates_path, 'locate_stats_res');
            end
        end
%         try 
%             load(save_stats_path);
%         catch
             %Feature statistics 
             [stats_features , dim_stat_types] = stats_main( features, algorithm_type, speech_labels, locate_stats_res);
              save(save_stats_path, 'stats_features', 'files_people_label', 'files_emotion_label', 'files_sentence_label');%save signal features
 %       end
        all_cross_feas = [all_cross_feas, stats_features];
 %       dim_stat_typess = [dim_stat_typess, dim_stat_types];
        stats_dims = [stats_dims, size(stats_features, 2)];
        %features, 
        if(issave(inx))
            save_stats_path1 = strcat(save_dir_main, feature_type{inx}, str1, '.mat');
            mdata = all_cross_feas;
            save(save_stats_path1, 'mdata', 'stats_dims', 'files_people_label', 'files_emotion_label', 'files_sentence_label');
            all_cross_feas = [];
        end
    end
%     catch
%     end
end
    %data_sel( all_cross_feas,  files_emotion_label);
end%�����������ڲ�ͬͳ����ʱʹ��
    end
end