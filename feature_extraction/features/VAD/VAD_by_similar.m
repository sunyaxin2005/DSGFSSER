function [ Distence_lfpc ] = VAD_by_similar( x, fs, frame_long, inc, alg_type, alg_prog )
%VAD_BY_SIMILAR Summary of this function goes here
%   Detailed explanation goes here
sun_lfpc_res = extract_features(x, fs, frame_long, inc, 'sun_lfpc', 20);
energys = sum(sun_lfpc_res, 2);
[B,IX] = sort(energys);
unvoice_len = size(sun_lfpc_res, 1)/10;
vector_noise = sum(sun_lfpc_res(IX(1 : unvoice_len), :));
vector_noise = vector_noise/unvoice_len;
for index = 1 : size(sun_lfpc_res, 1)
    Distence_lfpc = pdist2(vector_noise,sun_lfpc_res, 'euclidean');
end
[ Distence_lfpc] = wave_smooth(Distence_lfpc);
end

function [ wave_smooth_res] = wave_smooth(features)
[frame_num feature_dim] = size(features);
wave_smooth_res = [];
for frame_index = 1 : frame_num
   [c, l] = wavedec(features(frame_index, :), 2, 'db2');
   wave_smooth_res = [wave_smooth_res;c(1 : l(1))];
end
end