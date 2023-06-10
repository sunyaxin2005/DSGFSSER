function [ FeaIndexs ] = feature_selection_multi_main(  A, emotion_label, FeaNumCandi, sel_type, Para   )
%FEATURE_SELECTION_MULTI_MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
% Para.lambda = 0.1;
% Para.beta = 0.1;
options = [];
if ~isempty(emotion_label)
    options.gnd = emotion_label;
    options.nUseEigenfunction = length(unique(emotion_label));
else
    options.nUseEigenfunction = 10;
end
fea = A;
options.gnd = emotion_label;
 options.W =  Get_01Weight( fea, options.gnd );
 [Y, V] = get_low_dim1(A, options);
 
 switch(sel_type)
     case 'SPFS'
         [Ws] = OptimizationL21L22Common_M( A, Y, Para );
     case 'DSGFS'
         [ M ] = getLLEGRAPH( fea', 5 );
         [Ws] = OptimizationL21L22Common_manifod_M( A, Y, M, Para );
 end
 
%[Ws] = OptimizationL21L22Common_M( A, Y, Para );
FeaIndexs = cell(length(Ws), 1);
for ii = 1 : length(Ws)
    eigvec = Ws{ii};
    eigvec = max(abs(eigvec),[],2);
    [dump,idx] = sort(eigvec,'descend');
         
    FeaIndex = cell(1,length(FeaNumCandi));
    for i = 1:length(FeaNumCandi)   
        if FeaNumCandi(i) < length(idx)
            FeaIndex{i} = idx(1:FeaNumCandi(i));
        else
            FeaIndex{i} = idx(1:length(idx));
        end
    end
    FeaIndexs{ii} = FeaIndex;
end
end


function [Y, V] = get_low_dim1(fea, options)
if isfield(options,'W')
   W = options.W;
end
[V,D]=eig(W);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i) = V(:,index(i));
end
for i = 1:size(W,2)
    W(:,i) = W(:,i).*sqrt(sorted(i));
end
W = W(:, 1 : options.nUseEigenfunction);
Y = W;
V = sorted(1:options.nUseEigenfunction);
end