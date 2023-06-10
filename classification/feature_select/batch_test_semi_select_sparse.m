function [ output_args ] = batch_test_semi_select_sparse(wave_files_names, wave_files_path)
%BATCH_TEST_SEMI Summary of this function goes here
%   Detailed explanation goes here
file_num = length(wave_files_names);

save_dir = 'F:\dimension_reduction\';
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'locateLRGPSSDR', 'LRGAAFF', 'MMP' };
%ALg_types = {'LRGPSSDRCOS', 'LRGPSSDR', 'LRGPSSDR_SPARSE'};
%ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'PCA', 'LDA', 'locateLRGPSSDR'};
%ALg_types = {'LRGPSSDR', 'PCA', 'LDA'};
%ALg_types = {'mrmr', 'FS', 'LLFS', 'L21RFS', 'l1ls_featuresignM', 'l1ls_featuresignM1', 'l1ls_featuresign', 'SPFS', 'LLESPFS'};
%ALg_types = {'MCFS', 'mrmr', 'FS', 'LLFS', 'L21RFS', 'l1ls_featuresignM',  'l1ls_featuresignM1', 'SPFS', 'LLESPFS'};
%ALg_types = { 'FS', 'l1ls_featuresignM',  'l1ls_featuresignM1', 'SPFS', 'LLESPFS'};
%ALg_types = {'MCFS', 'mrmr', 'FS', 'LLFS','L21RFS','l1ls_featuresign', 'l1ls_featuresignM1', 'l1ls_featuresignM', 'SPFS', 'LLESPFS'};
ALg_types = {'MCFS','l1ls_featuresignM1',  'SPFS', 'LLESPFS'};
%ALg_types = { 'SemiFS', 'LSSFS', 'SLDFS', 'SLS', 'NaiveLS','MCFSM',  'MCFS'};%, 'condred','disr', 'jmi', 'cmim', 'PCA', 'LDA'};
lambda = [0.8];
beta = [0.5];
for ilambda = 1 : length(lambda)
    for ibeta = 1 : length(beta)
        Para.lambda = lambda(ilambda);
        Para.beta = beta(ibeta);
        fprintf('%3f_%3f\n', Para.lambda, Para.beta);
for ii = 1 : length(wave_files_names)
    fprintf('%s\n', wave_files_names(ii).name);
    for jj = 1 : length(ALg_types)
    %for jj = 1
        fprintf('%s\n', ALg_types{jj});
        %fprintf('%d\n', jj);
        
 %       try
           file_path = strcat(wave_files_path, wave_files_names(ii).name);
           X = [];
           load(file_path);
%            try 
%                size(X);
%            catch
               X = fea;
               Y = gnd;
 %          end
%             X = x';
%             Y = t';
           D = size(X, 2);
           if D > 200
               kk = 10:20:400;
           else
               step = floor(D/15);
               kk = 10 : step : D *2/3;
           end
            [ prs, pstds] = test_feature_selection_sparse( X,Y, ALg_types{jj}, kk, Para);
            
            
           %[ prs] = test_feature_selection_sparse( fea, gnd, ALg_types{jj}, kk);
          %[ prs] = test_feature_selection( fea, gnd, jj, kk);
           %[ prs] = test_semi( fea, gnd, ALg_types{jj}, kk);
           %prs = prs(1:1,:);
           prs = [prs, pstds];
           for aa = 1 : 5;%size(prs, 1)
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
end

end

