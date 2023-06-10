function [ output_args ] = batch_test_semi_select(wave_files_names, wave_files_path)
%BATCH_TEST_SEMI Summary of this function goes here
%   Detailed explanation goes here
file_num = length(wave_files_names);

save_dir = 'F:\dimension_reduction\';
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'locateLRGPSSDR', 'LRGAAFF', 'MMP' };
%ALg_types = {'LRGPSSDRCOS', 'LRGPSSDR', 'LRGPSSDR_SPARSE'};
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'PCA', 'LDA', 'locateLRGPSSDR'};
%ALg_types = {'LRGPSSDR', 'PCA', 'LDA'};
ALg_types = {'MCFSM'};
%ALg_types = { 'SemiFS', 'LSSFS', 'SLDFS', 'SLS', 'NaiveLS','MCFSM',  'MCFS'};%, 'condred','disr', 'jmi', 'cmim', 'PCA', 'LDA'};
for ii = 1 : length(wave_files_names)
    fprintf('%s\n', wave_files_names(ii).name);
    for jj = 1 : length(ALg_types)
    %for jj = 1
        fprintf('%s', ALg_types{jj});
        %fprintf('%d\n', jj);
        kk = 100:40:400;
        %try
           file_path = strcat(wave_files_path, wave_files_names(ii).name);
           load(file_path);
           [ prs] = test_feature_selection( fea, gnd, ALg_types{jj}, kk);
          %[ prs] = test_feature_selection( fea, gnd, jj, kk);
           %[ prs] = test_semi( fea, gnd, ALg_types{jj}, kk);
           %prs = prs(1:1,:);
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

