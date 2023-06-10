function [ output_args ] = batch_test_semi(wave_files_names, wave_files_path)
%BATCH_TEST_SEMI Summary of this function goes here
%   Detailed explanation goes here
file_num = length(wave_files_names);

save_dir = 'F:\dimension_reduction\';
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'locateLRGPSSDR', 'LRGAAFF', 'MMP' };
%ALg_types = {'LRGPSSDRCOS', 'LRGPSSDR', 'LRGPSSDR_SPARSE'};
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'PCA', 'LDA', 'locateLRGPSSDR'};
%ALg_types = {'LRGPSSDR',  'PCA', 'LDA', 'SDA', 'LSDA',  'NPE', 'LPP' };
%ALg_types = {'OPSRDR', 'DSPP_train', 'SPP_train', 'LRGPSSDR'};
ALg_types = { 'LDA','NPE', 'LPP' };
%ALg_types = {'KECAMCFS', 'MCFS'};
%ALg_types = {'MCFS', 'PCA', 'LDA', 'mrmr'};
global file_name; 
for ii = 1 : length(wave_files_names)
    fprintf('%s\n', wave_files_names(ii).name);
    for jj = 1 : length(ALg_types)
        fprintf('%s\n', ALg_types{jj});
        kk = 10:20:200;
        %kk = 100;
 %    try
       temp = wave_files_names(ii).name;
       file_name = strcat('F:\稀疏表示实验结果\img\', temp(1:3));
       %file_name = file_name(1:length(file_name);
           file_path = strcat(wave_files_path, wave_files_names(ii).name);
           load(file_path);
           %[ prs] = test_feature_selection( fea, gnd, ALg_types{jj}, kk);
           %[ prs, pstds] = locate_dimension_reduction( fea, gnd, ALg_types{jj}, kk);
            X = fea;
            Y = gnd;
%             X = x';
%             Y = t';
           [ prs, pstds] = test_semi( X, Y, ALg_types{jj}, kk);
           %prs = prs(1:1,:);
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
% %         
        save_path = strcat(save_dir, ALg_types{jj}, wave_files_names(ii).name);
        save(save_path, 'prs');
    end
end

end

