function [ M ] = get_sub_space_G_graph( fea, K )
%GET_SUB_SPACE_G_GRAPH Summary of this function goes here
%   Detailed explanation goes here
[N, D] = size(fea);
rand_times = 1;
M = 0;
for ii = 1 : rand_times
    s = sampling([1:D], floor(D/2));
    sub_fea = fea(:, s);
%     distance = pdist2(sub_fea, sub_fea);
%     [sorted,index] = sort(distance);
%     sita = sorted( K+1, :);
%     sita = sita' *  sita;
%     W = exp(-(distance.^2)./sita);
    [W1] = LLEM(sub_fea',K);
    M = M + W1;
end
M = M./rand_times;
end

