function [ output_args ] = batch_DSGFSSER( wave_files_names, wave_files_path )
%BATCH_HUWSF Summary of this function goes here
%   Detailed explanation goes here
%using example: cd F:\第一篇论文\berlin; res = dir('*.mat');batch_HuWSF(res, 'F:\第一篇论文\berlin\');
ALg_types = { 'DSGFS'}; %if you want to see the results of other feature selection, you can add {'condred', 'relief', 'mrmr', 'fcbf', 'disr', 'cmi', 'cmim', 'jmi', 'SFS'} into brace, for example:ALg_types = { 'MCFS', 'mrmr'}; 
lambda = [0.8];
beta = [0.5];
Para.lambda = lambda;
Para.beta = beta;
for ii = 1  : length(wave_files_names) 
    fprintf('%s\n', wave_files_names(ii).name);
    for jj = 1 : length(ALg_types)
    %for jj = 1
         fprintf('%s\n', ALg_types{jj});
        %fprintf('%d\n', jj);
        
       % try
           file_path = strcat(wave_files_path, wave_files_names(ii).name);
%           file_path1 = strcat(wave_files_path, wave_files_names(7).name);
           %file_path1 = 'F:\第一篇论文\datas\bSAVEE_OTHER';
           %file_path1 = 'F:\第一篇论文\datas\bSAVEE_OTHER';
           %file_path1 = 'F:\第一篇论文\berlin\berlinOTHER';
%           load(file_path1);
           %X = [];
           load(file_path);
           %X = all_cross_feas;
           
            people_label =  files_people_label;
            X = mdata;
            Y = files_emotion_label;
%             X = fea;
%             Y  = gnd;

           D = size(X, 2);
           if ii <= 5
               high_threhold = 400;
               if D > high_threhold
                   high_threhold = 400;
               else
                   high_threhold = D;
               end
                kk = 10 : 40 : high_threhold;
               %is_pre_sel = 0;
               pre_sel_threhold = 0;
           else
               high_threhold = 600;
               kk = 10:20:high_threhold;
               if length(unique(people_label)) == 10
                    pre_sel_threhold = 0.10;% pre_sel_threhold = 0.10, for SD on EmoDB 
               else 
                    pre_sel_threhold = 0.05;%pre_sel_threhold = 0.05, for SD on SAVEE and CASIA 
               end   
           end
pre_sel_threhold = 0.05;
            [ prs, pstds] = test_DSGFSSER( X,Y, people_label, ALg_types{jj}, kk, Para, 'PI', pre_sel_threhold);
            
           prs = [prs, pstds];
           for aa = 1 : size(prs, 1)
               fprintf('alg%d:', aa);
               for bb = 1 : size(prs, 2)
                   fprintf(' %f', prs(aa, bb));
               end
               fprintf('\n', aa);
           end
%         catch
%          prs = [];
%         end
            
%         save_path = strcat(save_dir, ALg_types{jj}, wave_files_names(ii).name);
%         save(save_path, 'prs');
    end

end
end

function [ prs, pstds] = test_DSGFSSER(X, L, people_label, alg_type, sel_dims, Para,  cross_type, pre_sel_threhold)
if strcmp(cross_type, 'PI') % speaker independent 
    random_times= 1;
    cross_type = 'LOSO';
    pre_sel_threhold = 0.05; % because for 
else% speaker dependent 
    cross_type = 'five_cross';
    random_times = 10;
end
%[ X ] = data_maxminnorm( X );
[ X ] = sun_norm_by_mean_std( X );
X(isnan(X)) = 0;
X(isinf(X)) = 0;
X = normalizated_by_people_label(X, people_label);
X(isnan(X)) = 0;
X(isinf(X)) = 0;
classifier_type = {'SVM'};
%classifier_type = {'SVM'};
pall = zeros(length(classifier_type) * 4 , random_times, length(sel_dims) + 1);
train_inxs_num = 0;
for ii = 1 : random_times
    train_inxs = get_train_inxs(cross_type, L, people_label); 
    train_inxs_num = length(train_inxs);
    for kk = 1 : length(train_inxs)
        train_inx = train_inxs{kk};
        XLTrain = X(train_inx, :);
        LTrain = L(train_inx);
        XTest = X(~train_inx, :);
        LTest = L(~train_inx);
% %     
       if  pre_sel_threhold > 0
           FDRs = FDR(XLTrain, LTrain );
           XLTrain = XLTrain(:, FDRs > pre_sel_threhold);
           XTest = XTest(:, FDRs > pre_sel_threhold);
       end
        
        [FeaIndexs ] = feature_selection_multi_main( XLTrain, LTrain, sel_dims, alg_type, Para);
        for jj = 1 : length(sel_dims)
            classesss = [];
            for kkkk = 1 : length(FeaIndexs)
                FeaIndex = FeaIndexs{kkkk};
                XLTrain2 = XLTrain(:, FeaIndex{jj});
                XTest2 = XTest(:, FeaIndex{jj});      
                [aa,  pc, classes] = classifier_main( XLTrain2, XTest2, LTrain, LTest, classifier_type{1});
                classesss = [classesss;classes(:)'];
            end
            [ vote_res, max_votes_num ] = vote_by_result( classesss );
            classes = vote_res(:);

            qq = 1;     
            %[compute_accuracy( classes, LTest ), compute_weight_recall( classes, LTest ), compute_unweight_recall( classes, LTest ), computer_percision( classes, LTest )]
            pall((qq - 1) * 4 + 1, ii, jj) = pall((qq - 1) * 4 + 1, ii, jj) + compute_accuracy( classes, LTest );%计算带权召回率，compute weighted recall
            pall((qq - 1) * 4 + 2, ii, jj) = pall((qq - 1) * 4 + 2, ii, jj) + compute_weight_recall( classes, LTest );%计算不带权召回率，compute unweighted recall
            pall((qq - 1) * 4 + 3, ii, jj) = pall((qq - 1) * 4 + 3, ii, jj) + compute_unweight_recall( classes, LTest );%计算不带权召回率，compute unweighted recall
            pall((qq - 1) * 4 + 4, ii, jj) = pall((qq - 1) * 4 + 4, ii, jj) + computer_percision( classes, LTest );%计算精度，computer percision
            
        end 
    end
   
end
 pall = pall./ train_inxs_num;

%计算交叉验证的方差, 均值， compute the mean and std for cross 
pii = zeros(random_times, length(sel_dims) + 1);
prs = [];
pstds = [];
for ii = 1 : size(pall, 1)
    pii(:, :) = pall(ii, :, :);
    p1s = std(pii, 0, 1);
    p1 = mean(pii, 1);
    p1(length(sel_dims) + 1) = max(p1);
    [~, inxxx] = max(p1);
    p1s(length(sel_dims) + 1) = p1s(inxxx);
    prs = [prs;p1];
    pstds = [pstds;p1s];
end
end

function X = normalizated_by_people_label(X, people_label)
UP = unique(people_label);
for inx = 1 :length(UP)
%     XII = X(people_label == UP(inx), :);
%     p_II = people_label(people_label == UP(inx), :);
%     inxs = selected_normalized(p_II, 21);
%     
%     mean_p = mean(XII(inxs, :));
%     std_p = std(XII(inxs, :));
%     mean_p1 = repmat(mean_p, size(XII, 1), 1);
%     std_p1 = repmat(std_p, size(XII, 1), 1);
%     XII = (XII - mean_p1)./1;
%     
    %X1(people_label == UP(inx), :) = XII;
    mean_p = mean(X(people_label == UP(inx), :));
    std_p = std(X(people_label == UP(inx), :));
    mean_p1 = repmat(mean_p, length(people_label(people_label == UP(inx))), 1);
    std_p1 = repmat(std_p, length(people_label(people_label == UP(inx))), 1);
    std_p1(std_p1<1e-4) = 1;
    X(people_label == UP(inx), :) = (X(people_label == UP(inx), :) - mean_p1)./1;
end
end

function [ vote_res, max_votes_num ] = vote_by_result( reco_result )
%VOTE_BY_RESULT Summary of this function goes here
%   Detailed explanation goes here
[ticket_num, sample_num] = size(reco_result);
%vote_res = [];
for samle_inx = 1 : sample_num
    fea = reco_result(:, samle_inx);
    u = unique(fea);
    max_num = 0;
    max_label = 0;
    for inx = 1 : length(u)
        this_num = length (fea(fea == u(inx)));
        if this_num > max_num
            max_num = this_num;
            max_label = u(inx);
        end
    end
    vote_res(samle_inx) = max_label;
    max_votes_num(samle_inx) = max_num;
end

end

