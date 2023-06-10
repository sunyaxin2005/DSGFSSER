function [ features ] = extract_features( x, fs, frame_long, inc, alg_type, alg_prog )
%EXTRACT_FEATURES Summary of this function goes here
%   Detailed explanation goes here
% x,单通道声音数据（sample_num * 1）
% fs,采样率
% frame_long，帧长度
% inc, 帧之间重叠的比例
% alg_type, 算法类型, 目前有：'short_energy', 'short_zerocross', 'short_correct',
% 'short_fft', 'ausees_erb', 'LPC', 'LPCC','LPCC_ACW', 'LPCC_PLF' 'MFCC','WMFCC' 'sun_mfcc', 'sun_lfpc', 'ZCPA', 'fractal_dim',
% 'ausees_bark', 'IMF', 'mel_bank_filter', 'HNR', 'TEO', 'MSF', 'PLP',
% 'harmonic', 'LFPC', 'MFCC_LOW', 'LP_MFCC', 'ZCMT', 'pitch',
% 'spectrl_energy', 'jitter', 'shimmer', 'mspectrum', 'rspectrum',
% 'sun_rspectrum', 'sun1_rspectrum', 'msf', 'lrrl', 'grad_angle', 'lbp',
% 'hu','dmqi', 'ver_grad', 'grad_in', 'hu1', 'entropy', 'lrrl_dct',
% 'grad_angle_dct', 'lbp_dct', 'hu_dct', 'dmqi_dct', 'ver_grad_dct',
% 'grad_in_dct', 'hu1_dct', 'entropy_dct','hu_ver_grad_dct'
%'unframe_mfcc'
% alg_prog， 算法参数
%以下需要使用到mir包
switch alg_type
    case 'short_energy'
        %short_energy_res, frame_num * 1
        features = short_energy(x, fs, frame_long, inc, alg_prog);
    case 'short_zerocross'
        %zero_cross_res, frame_num * 1
        features = short_zerocross(x, fs, frame_long, inc, alg_prog);
    case 'short_correct'
        %short_energy_res, frame_num * frame_long/2 
        features = short_correct(x, fs, frame_long, inc, alg_prog);
    case 'short_fft'
        %short_fft_res, frame_num * frame_long向下最小的2的N次方+1
        features = short_fft(x, fs, frame_long, inc, alg_prog);
    case 'ausees_bark'
        %ausees_erb_res: frame_num * 17
        features = ausees_bark(x, fs, frame_long, inc, alg_prog);
    case 'ausees_erb'
         %ausees_erb_res: frame_num * 17
        features = ausees_erb(x, fs, frame_long, inc, alg_prog);
    case 'LPC'
        %frame_num * p
        features = LPC(x, fs, frame_long, inc, 'LPC');
%         temp_features = features;
%         features = [];
%         for index = 1 : size(temp_features, 1)
%             temp = temp_features(index, :);
%             features = [features;dct(temp)];
%         end
    case 'LPCC'
        %frame_num * p
        features = LPC(x, fs, frame_long, inc, 'LPCC');
    case 'LPCC_ACW'
        %frame_num * p
        features = LPC(x, fs, frame_long, inc, 'LPCC_ACW');
    case 'LPCC_PLF'
        %frame_num * p
        features = LPC(x, fs, frame_long, inc, 'LPCC_PLF');     
    case 'MFCC'
        %frame_num * 13
        features = mfcc(x, fs, frame_long, inc, alg_prog);
    case 'WMFCC'
        features = WMFCC(x, fs, frame_long, inc, alg_prog);
   case 'unframe_mfcc'
        %frame_num * 13
        features = unframe_mfcc(x, fs, frame_long, inc, alg_prog);       
   case 'sun_mfcc'
        %frame_num * 12
        features = sun_mfcc(x, fs, frame_long, inc, alg_prog);
  case 'sun_lfpc'
        %frame_num * 12
        features = sun_lfpc(x, fs, frame_long, inc, alg_prog);
%'lrrl', 'grad_angle', 'lbp',
% 'hu','dmqi', 'ver_grad', 'grad_in', 'hu1', 'entropy', 'lrrl_dct',
% 'grad_angle_dct', 'lbp_dct', 'hu_dct', 'dmqi_dct', 'ver_grad_dct',
% 'grad_in_dct', 'hu1_dct', 'entropy_dct'
    case 'lrrl'
        features = f_lrrl(x, fs, frame_long, inc, alg_prog);
     case 'grad_angle'
        features = f_grad_angle(x, fs, frame_long, inc, alg_prog); 
    case 'lbp'
        features = f_lbp(x, fs, frame_long, inc, alg_prog); 
    case 'Q'
        features = f_Q(x, fs, frame_long, inc, alg_prog); 
     case 'dmqi'
        features = f_dmqi(x, fs, frame_long, inc, alg_prog);
     case 'ver_grad'
        features = f_ver_grad(x, fs, frame_long, inc, alg_prog); 
    case 'grad_in'
        features = f_grad_in(x, fs, frame_long, inc, alg_prog); 
    case 'hu1'
        features = f_hu1(x, fs, frame_long, inc, alg_prog);     
    case 'entropy'
        features = f_entropy(x, fs, frame_long, inc, alg_prog);    
    case 'mean'
        features = f_mean(x, fs, frame_long, inc, alg_prog);          
    case 'lrrl_dct'
        features = f_lrrl_dct(x, fs, frame_long, inc, alg_prog);
     case 'grad_angle_dct'
        features = f_grad_angle_dct(x, fs, frame_long, inc, alg_prog); 
    case 'lbp_dct'
        features = f_lbp_dct(x, fs, frame_long, inc, alg_prog); 
    case 'TQ'
        features = f_TQ(x, fs, frame_long, inc, alg_prog); 
     case 'dmqi_dct'
        features = f_dmqi_dct(x, fs, frame_long, inc, alg_prog);
     case 'ver_grad_dct'
        features = f_ver_grad_dct(x, fs, frame_long, inc, alg_prog); 
    case 'grad_in_dct'
        features = f_grad_in_dct(x, fs, frame_long, inc, alg_prog); 
    case 'hu_ver_grad_dct'
        features = f_hu_ver_grad_dct(x, fs, frame_long, inc, alg_prog);         
    case 'hu1_dct'
        features = f_hu1_dct(x, fs, frame_long, inc, alg_prog);     
    case 'entropy_dct'
        features = f_entropy_dct(x, fs, frame_long, inc, alg_prog);   
    case 'mean_dct'
        features = f_mean_dct(x, fs, frame_long, inc, alg_prog);            
    case 'ZCPA'
        %frame_num * 16
        features = ZCPA( x,  fs, frame_long, inc, 10);
    case 'fractal_dim'
        %frame_num * 1
        features = fractal_dim(x,  fs, frame_long, inc, alg_prog);
    case 'IMF'
        %sample_num * imf_coefs
        features = imf(x,  fs, frame_long, inc, alg_prog);
    case 'mel_bank_filter'
        %frame_num * 37
        features = mel_bank_filter(x,  fs, frame_long, inc, alg_prog);
    case 'HNR'
        %frame * 1
    case 'TEO'
        %sample_num * 1
        features = TEO(x,  fs, frame_long, inc, alg_prog);
    case 'MSF'
    case 'PLP'
        %frame * p
    case 'harmonic'
        %51 * 14, fft频率数 * 滤波器数,全局参数
        features = harmonic(x,  fs, frame_long, inc, alg_prog);
    case 'LFPC'
        %frame_num * 12, 12维滤波器数
        features = LFPC(x,  fs, frame_long, inc, alg_prog);
    case 'MFCC_LOW'
        
    case 'LP_MFCC'
        
    case 'ZCMT'
        %frame_num * 16
        features = ZCMT( x,  fs, frame_long, inc, 10);
    case 'pitch'
        %frame_num * 1
        features = pitch( x,  fs, frame_long, inc, alg_prog);
    case 'spectrl_energy'
        %frame_num * 15, 15为滤波器数
        features = spectrl_energy( x,  fs, frame_long, inc, alg_prog);
    case 'jitter'
        %frame * 1
        features = jitter( x,  fs, frame_long, inc, alg_prog);
    case 'shimmer'
        %frame * 1
        features = shimmer( x,  fs, frame_long, inc, alg_prog);
    case 'mspectrum'
        features = spectrum( x,  fs, frame_long, inc, alg_prog);
    case 'rspectrum'
        features = rspectrum( x,  fs, frame_long, inc, alg_prog);
    case 'sun_rspectrum'
        features = sun_rspectrum( x,  fs, frame_long, inc, alg_prog);
    case 'sun1_rspectrum'
        features = sun1_rspectrum( x,  fs, frame_long, inc, alg_prog);        
    case 'msf'
        features = msf(x, fs, frame_long, inc, alg_prog );
end
end

function [ short_energy_res ] = short_energy( x, fs, frame_long, inc, alg_prog )
%short_energy_res, 1 * frame_num
%x = awgn(x, 20);
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
rms = mirrms(f);
short_energy_res = mirgetdata(rms);
short_energy_res = short_energy_res';
end

function [ zero_cross_res ] = short_zerocross( x, fs, frame_long, inc, alg_prog )
%zero_cross_res, 1 * frame_num
x = awgn(x, 35);
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
zero_cross = mirzerocross(f);
zero_cross_res = mirgetdata(zero_cross);
zero_cross_res = zero_cross_res';
end

function [ short_correct_res ] = short_correct( x, fs, frame_long, inc, alg_prog )
%short_energy_res, frame_long/2 * frame_num
x = awgn(x, 35);
%figure, plot(x);
% fft_x = fft(x);
% fft_x(length(fft_x)/8 : length(fft_x) * 7/8) = 0;
% ifft_x = ifft(fft_x);
% x = real(ifft_x);
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
correct_res = mirautocor(f);
short_correct_res = mirgetdata(correct_res);
short_correct_res = short_correct_res';
end

function [ short_fft_res ] = short_fft( x, fs, frame_long, inc, alg_prog )
%short_fft_res, frame_long向下最小的2的N次方+1 * frame_num
%为什么如此见mirspectrum的分析
% a = miraudio(x, fs);
% f = mirframe(a, frame_long, 'sp', inc,  'sp');
% f = mirgetdata(f);
% [ sample_num frame_num] = size(f);
% short_fft_res = [];
% for frame_index = 1 : frame_num
%     %figure,plot(f(:, frame_index));
%     %lpc_trans_res = lpc_trans(f(:, frame_index)');
%     %xxx = f(:, frame_index);
%     [ window_res ] = hanning_window( f(:, frame_index) );
%     %figure,plot(f(:, frame_index)' - lpc_trans_res);    
%     HX = fft(window_res');
%     short_fft_res = [short_fft_res ;abs(HX(1 : length(HX)/2))];
%     %figure,plot(abs(HX(1 : length(HX)/2)));
% end
% 
% % 
% %  fft_res = mirspectrum(f);
% %  short_fft_res = mirgetdata(fft_res);
%  short_fft_res = log10(short_fft_res);
% x = awgn(x, 40);
 [T,F,B] = spgrambw(x(:, 1),fs,'pjcwm', 60);
 short_fft_res =log10(B);
% meanb = min(min(short_fft_res));
% short_fft_res = short_fft_res - meanb;
% [frame_num frequency_num] = size(B);
% for frame_index = 1 : frame_num
%     short_fft_res(frame_index, :) = TEO(short_fft_res(frame_index, :)')';
% end

end

function [ lpc_trans_res ] = lpc_trans(feature)
[frame_num feature_dim] = size(feature);
lpc_trans_res = [];
p_coef_num = 4;
for frame_index = 1 : frame_num
    ar = lpc(feature(frame_index, :), p_coef_num);
    lpc_value = feature(frame_index, :);
    for frame_dim_index = p_coef_num + 1 : feature_dim
        lpc_value(1, frame_dim_index) = 0;
        for p_coef_index = 1 : p_coef_num 
            lpc_value(1, frame_dim_index) = lpc_value(1, frame_dim_index) - feature(frame_index, frame_dim_index - p_coef_index) * ar(1, p_coef_index + 1);
        end
    end
    lpc_trans_res = [lpc_trans_res;lpc_value];
end
end

function [ ausees_erb_res ] = ausees_bark( x, fs, frame_long, inc, alg_prog )
%ausees_erb_res: frame_num * 17
%frame_long: 256, 帧移：128
%第一步：预加重（高通滤波器），分帧，加窗（hanning）
%第二步：FFT变换，得到频谱上的能量分布
%第三步：对每一帧信号计算其幅度的平方，然后取对数，得到对数能量谱
%第四步：域值判断，将0db作为域值
%第五步：频带划分，按照bark频带划分，只选前17个
%第六步：对于每个子带，分别计算其谱包络下的面积
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[ sample_num frame_num] = size(f);

ausees_erb_res = [];
for k = 1 : frame_num
    window_frame = hanning_window(f(:,  k));%第一步：预加重（高通滤波器），分帧，加窗（hanning）
    fft_res = fft(window_frame,sample_num,1);%第二步：FFT变换，得到频谱上的能量分布
%     ample_res = abs(fft_res);%第三步：对每一帧信号计算其幅度
%     for sample_index = 1 : sample_num
%         if log10(ample_res(sample_index) * ample_res(sample_index)) <= 0 %第四步：域值判断，将0db作为域值
%             fft_res(sample_index) = 0;
%         end
%     end
    y = ifft(fft_res, sample_num);
    frame_a = miraudio(y, fs);
    frame_bank = mirfilterbank(frame_a, 'Bark');%19
    enve_res = mirenvelope(frame_bank, 'PostDecim',1);%求包络
    enve_res = mirgetdata(enve_res);
    [frame_sample_num frame_frame_num filter_num] = size(enve_res);
    areas = zeros(1, 19);
    
    for ps_index = 1 : 19 %对于每个子带，分别计算其谱包络下的面积
        if ps_index > filter_num %mir包进行了优化，有可能有些的信号没有那么多带通滤波结果
            continue;
        end
        for frame_sample_index = 1: frame_sample_num
            areas(1, ps_index) = areas(1, ps_index) + enve_res(frame_sample_index, 1, ps_index);
        end
    end
    
    ausees_erb_res = [ausees_erb_res;areas];
end
    
end

function [ ausees_erb_res ] = ausees_erb( x, fs, frame_long, inc, alg_prog )
%frame_long: 256, 帧移：128
%第一步：预加重（高通滤波器），分帧，加窗（hanning）
%第二步：FFT变换，得到频谱上的能量分布
%第三步：对每一帧信号计算其幅度的平方，然后取对数，得到对数能量谱
%第四步：域值判断，将0db作为域值
%第五步：频带划分，按照ERB频带划分，只选前17个, 与ausees_bark仅仅这个地方不一样
%第六步：对于每个子带，分别计算其谱包络下的面积
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[ sample_num frame_num] = size(f);

ausees_erb_res = [];
for k = 1 : frame_num
    window_frame = hanning_window(f(:,  k));%第一步：预加重（高通滤波器），分帧，加窗（hanning）
    fft_res = fft(window_frame,sample_num,1);%第二步：FFT变换，得到频谱上的能量分布
    ample_res = abs(fft_res);%第三步：对每一帧信号计算其幅度
    for sample_index = 1 : sample_num
        if log10(ample_res(sample_index) * ample_res(sample_index)) <= 0 %第四步：域值判断，将0db作为域值
            fft_res(sample_index) = 0;
        end
    end
    y = ifft(fft_res, sample_num);
    frame_a = miraudio(y, fs);
    frame_bank = mirfilterbank(frame_a, 'Gammatone', 'NbChannels', 17);%17
    enve_res = mirenvelope(frame_bank, 'PostDecim',1);%求包络
    enve_res = mirgetdata(enve_res);
    [frame_sample_num frame_frame_num filter_num] = size(enve_res);
    areas = zeros(1, 17);
    
    for ps_index = 1 : 17 %对于每个子带，分别计算其谱包络下的面积
        if ps_index > filter_num %mir包进行了优化，有可能有些的信号没有那么多带通滤波结果
            continue;
        end
        for frame_sample_index = 1: frame_sample_num
            areas(1, ps_index) = areas(1, ps_index) + enve_res(frame_sample_index, 1, ps_index);
        end
    end
    
    ausees_erb_res = [ausees_erb_res; areas];
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lpcc及其变形系数求解区
function [ lpc_res] = LPC( x, fs, frame_long, inc, alg_prog)
%LPCC Summary of this function goes here
%   Detailed explanation goes here
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[ sample_num frame_num] = size(f);
lpc_res = [];
for frame_index = 1 : 1:frame_num
    ar = lpcauto(f(:, frame_index));
    switch alg_prog
        case 'LPC'
           [row, p] = size(ar);
           ar = ar(1, 2 : p);
           lpc_res = [lpc_res, ar(:)];
        case 'LPCC'
           c = LPCC( ar);
           lpc_res = [lpc_res, c(:)];
        case 'LPCC_ACW'
            cacw = ACW( ar);
            lpc_res =  [lpc_res, cacw(:)];
        case 'LPCC_PLF'
            cplf = PLF( ar);
            lpc_res =  [lpc_res, cplf(:)];    
    end
end
lpc_res = lpc_res';
end

function [ c ] = LPCC( ar)
[row, p] = size(ar);
ar = ar(1, 2 : p);%因为计算LPC的库在ar第一个位置上人为的添加了1，所以去掉
p = p - 1;
c = zeros(1, p);
c(1) = ar(1);
for n = 2 : p
    x = 0;
    for k = 1 : n - 1
        x = x + (1 - k/n) * ar(k) * c(n - k);
    end
    c(n) = ar(n) - x;
end

end

function [ cacw ] = ACW( ar)
clp = LPCC( ar);
[row, p] = size(ar);
ar = ar(1, 2 : p);%因为计算LPC的库在ar第一个位置上人为的添加了1，所以去掉
p = p - 1;
cnn = zeros(1, p);
cacw = zeros(1, p);

b = (p - 1) * ar(1)/p;
cnn(1) = b;
cacw(1) = log(p);
for n = 2 : p
    x = 0;
    for k = 1 : n - 1
        b = (p - k) * ar(k) / p;
        x = x + (1 - k/n) * b * cnn(n - k);
    end
    cnn(n) = (p - n) * ar(n) / p - x;
    cacw(n) = clp(n) -  cnn(n);
end

end

function [ cplf ] = PLF( ar)
clp = LPCC( ar);
[row, p] = size(clp);
cplf = zeros(1, p);
a = 1.0;
b = 0.9;
for n = 1 : p
    cplf(n) = clp(n) * (a^n - b^n);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mfcc_res ] = mfcc( x, fs, frame_long, inc, alg_prog )
%mfcc, frame_num * mfcc_coefs
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
mfcc_res = mirmfcc(f);
mfcc_res = mirgetdata(mfcc_res);
mfcc_res = mfcc_res';
mfcc_res = mfcc_res(:, 2:13);
end
 function [ mfcc_res ]  = WMFCC(x, fs, frame_long, inc, alg_prog)
 [ mfcc_res ] = wsf( x, fs, frame_long, inc, alg_prog );
 end
function [ mfcc_res ] = unframe_mfcc( x, fs, frame_long, inc, alg_prog )
%mfcc, frame_num * mfcc_coefs
a = miraudio(x, fs);
mfcc_res = mirmfcc(a);
mfcc_res = mirgetdata(mfcc_res);
end
function [ b ] = sun_mfcc( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = sun_lfpc( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
end

function [ b ] = f_lrrl( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = sqrt(lr_grad .^2 + rl_grad.^2);
end
function [ b ] = f_grad_angle( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = grad_angle;
end
function [ b ] = f_lbp( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
            mapping=getmapping(8, 'u2'); 
            [b,CLBP_MH]=CLBP(b,1,8,mapping,'i');
end
function [ b ] = f_Q( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = speech_hu( b, 5 );
b = b1;
end
function [ b ] = f_dmqi( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ b denoise] = dmqi( b );
end
function [ b ] = f_ver_grad( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = ver_grad;
end
function [ b ] = f_grad_in( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = grad_intenisty;
end
function [ b ] = f_hu1( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = img_hu( b );
b = b2;
end
function [ b ] = f_entropy( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
b = im_entropy(b);
end
function [ b ] = f_mean( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
b = im_mean(b);
end
function [ b ] = f_lrrl_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = sqrt(lr_grad .^2 + rl_grad.^2);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_grad_angle_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  , ~, lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = grad_angle;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_lbp_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
b = LBP_RAW(b);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_TQ( x, fs, frame_long, inc, alg_prog )
%x = awgn(x, 30);
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = speech_hu( b, 5 );
b = b1;
TQ = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    TQ = [TQ; dct_res(2 : 13)];
end
b = TQ;
end
function [ b ] = f_dmqi_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ b denoise] = dmqi( b );
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_ver_grad_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = ver_grad;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end

function [ b ] = f_hu_ver_grad_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = img_hu( b );
b = b1;
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = ver_grad;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end

function [ b ] = f_grad_in_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[ hor_grad  ver_grad lr_grad rl_grad grad_intenisty grad_angle] = direct_field(b); 
b = grad_intenisty;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end

function [ b ] = f_hu1_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = img_hu( b );
b = b2;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_entropy_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
b = im_entropy(b);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ b ] = f_mean_dct( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
b = im_mean(b);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: end)];
end
b = mfcc_res;
end
function [ y ] = ZCPA( X,  fs, frame_long, inc, filter_bank_num)
%ZCPA Summary of this function goes here
%  frame_num * 16
%因为需要采样率，所以X是matlab自带的wavread读取到
%ZCPA程序步骤
%第一步：分帧
%第二步：滤波
%第三步：对每个滤波结果，对每个帧计算短时过零率
%第四步：通过过零率得到箱号
%第五步：计算直方图
X = miraudio(X(:), fs);

% f = mirgetdata(f);
% [sample_num frame_num] = size(f);
frequency_bin = [0, 106, 223, 352, 495, 655, 829, 1022, 1236, 1473,...
    1734, 2024, 2344, 2689, 3089, 3522, 4000];

filter_res = mirfilterbank(X,  'NbChannels', filter_bank_num);
f = mirframe(filter_res, frame_long, 'sp', inc,  'sp');
%f
filter_res = mirgetdata(f );
[sample_num frame_num filter_bank_num] = size(filter_res);
y = zeros(frame_num, 16);
insss = [];
for filter_bank_index = 1 : filter_bank_num
    for frame_index = 1 : frame_num           
        [ peaks peaks_pos up_zerocross] = find_peak_between_zeros(filter_res(:, frame_index, filter_bank_index));
        zerocross = length(up_zerocross);
        zerocross = zerocross * fs/frame_long;
        if ~isempty(peaks)
            peaks_sum = sum(log10(peaks))/length(peaks);
        else
            peaks_sum = 0;
        end
        
        w = zerocross(1, 1);
        in = find(frequency_bin > w);
        [row col] = size(in);
        if(col > 0)
            bin_num = in(1, 1);
            bin_num = bin_num - 1;
            if bin_num < 17
                y(frame_index, bin_num) = y(frame_index, bin_num) + peaks_sum;
            end
        end
    end
end

end

function [ y ] = ZCMT( X,  fs, frame_long, inc, filter_bank_num)
%ZCPA Summary of this function goes here
%   Detailed explanation goes here
%因为需要采样率，所以X是matlab自带的wavread读取到
%ZCPA程序步骤
%第一步：分帧
%第二步：滤波
%第三步：对每个滤波结果，对每个帧计算短时过零率
%第四步：通过过零率得到箱号
%第五步：计算直方图
X = miraudio(X(:), fs);
% f = mirgetdata(f);
% [sample_num frame_num] = size(f);
frequency_bin = [0, 106, 223, 352, 495, 655, 829, 1022, 1236, 1473,...
    1734, 2024, 2344, 2689, 3089, 3522, 4000];

filter_res = mirfilterbank(X,  'NbChannels', filter_bank_num);
f = mirframe(filter_res, frame_long, 'sp', inc,  'sp');
filter_res = mirgetdata(f );
[sample_num frame_num filter_bank_num] = size(filter_res);
y = zeros(frame_num, 16);

for filter_bank_index = 1 : filter_bank_num
    for frame_index = 1 : frame_num
        [ peaks up_zerocross down_zerocross] = find_max_teo_between_zeros(filter_res(:, frame_index, filter_bank_index));%与ZCPT旧这个地方不一样，就变成了一篇SCI
        zerocross = length(up_zerocross);
        zerocross = zerocross * fs/frame_long;
        w = zerocross(1, 1);
        in = find(frequency_bin > w);
        [row col] = size(in);
        if(col > 0)
            bin_num = in(1, 1);
            bin_num = bin_num - 1;
            if bin_num < 17
                peaks_sum = sum(peaks);
                y(frame_index, bin_num) = y(frame_index, bin_num) + peaks_sum;
            end
        end
    end
end

end

function [ hw_res ] = fractal_dim( x, fs, frame_long, inc, alg_prog )
%hw_res, frame_num * 1
%
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[sample_num frame_num] = size(f);
hw_res = zeros(frame_num, 1);
for frame_index = 1 : frame_num
     hw_res(frame_index, 1) = hw_boxdim(f(:,  frame_index)', 10);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imf相关系数提取区域
function [ imf_res ] = imf( x, fs, frame_long, inc, alg_prog )
%imf_res, sample_num * imf_coefS
%
imf = emd(x);
M = length(imf);
[sample_num frame_num] = size( x); 
imf_res = zeros(sample_num * M);

for k = 1 : M
    imf_res = [imf_res, imf{k}];
%     [Hkin hil_res] = hilbert_envelope(imf{k});
%     Q = atan(hil_res/imf{k});
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mel_bank_filter_res ] = mel_bank_filter( x, fs, frame_long, inc, alg_prog )
%hw_res, sample_num * imf_coefS
%mel滤波，然后分帧，最后求每一帧的能量
a = miraudio(x, fs);
filter_res = mirfilterbank(a, 'Mel');
f = mirframe(filter_res, frame_long, 'sp', inc,  'sp');
rms = mirrms(f);
rms = mirgetdata(rms);
mel_bank_filter_res = [];
[sample_num frame_num filter_num] = size(rms);
temp = zeros(1, filter_num);
for frame_index = 1 : frame_num
    temp(1, :)  = rms(1, frame_index, :);
    mel_bank_filter_res = [mel_bank_filter_res;temp(1, :)];
end

end

function [ teo_res ] = TEO( x, fs, frame_long, inc, alg_prog )
%hw_res, sample_num * 1
%mel滤波，然后分帧，最后求每一帧的能量
[sample_num channel_num] = size(x);

teo_res = zeros(sample_num, channel_num);
for channel_index = 1: channel_num
    for sample_index = 2 : sample_num - 1
        teo_res(sample_index, channel_index) = x(sample_index, channel_index)^2 - x(sample_index - 1, channel_index) * x(sample_index + 1, channel_index);
    end
end

end

function [ harmonic_res ] = harmonic(  x, fs, frame_long, inc, alg_prog  )
%HARMONIC Summary of this function goes here
%   harmonic_res, 51 * 14, fft频率数 * 滤波器数
a = miraudio(x, fs);
%第一步分帧
% x
% f = mirframe(x, frame_long, 'sp', frame_long,  'sp');
% df = mirgetdata(f);
% %第二步求没一帧的基音频率
% p = mirpitch(f, 'AutocorSpectrum');
% dp = mirgetdata(p);
% %第三步Time varying subband filter
% [ sample_num frame_num] = size(df);
% total_filter_res = [];
% zeroframe = zeros(sample_num, 14);
% for frame_index = 1 : frame_num
%     if(~isnan(dp(1, frame_index)) && dp(1, frame_index) > 20)  
%         fbank = (dp(1, frame_index)/2 : dp(1, frame_index) : dp(1, frame_index) * 15);
%         fa = miraudio(df(:, frame_index));
%         filter_res = mirfilterbank(fa, 'Manual', fbank);
%         filter_res = mirgetdata(filter_res);
%         dfilter_res = [];
%         [p1 q1 r1] = size(filter_res);
%         if r1 ~= 14
%             continue;
%         end
%         for filter_index = 1 : 14
%             dfilter_res = [dfilter_res, filter_res(:, 1, filter_index)];
%         end   
%         total_filter_res = [total_filter_res;dfilter_res];
%     else
%         total_filter_res = [total_filter_res;zeroframe];
%     end
% end
total_filter_res = mirfilterbank(a, 'NbChannels', 14);
total_filter_res = mirgetdata(total_filter_res);
 % 第四步得到每一滤波后信号的上包络, 并求包络的快速傅立叶变换
 harmonic_res = [];
 for filter_index = 1 : 14
    a = miraudio(total_filter_res(:, filter_index), fs);
   %e = mirenvelope(a,'Gauss',2, 'Smooth', 2, 'PostDecim',1);
    a = mirgetdata(a);
    e = hilbert_envelope(a);
    fft_res = rfft(e, 100, 1);
    fft_res = abs(fft_res(:, 1));
    harmonic_res = [harmonic_res, fft_res];
 end

end

function [ set ] = LFPC(  x, fs, frame_long, inc, alg_prog  )
%LFPC Log Frequency Power Coefficients
%   Detailed explanation goes here
a = miraudio(x, fs);
filter_bank = [200, 250, 315, 400, 500, 630, 800, 1000, 1260, 1600, 2000, 2520, 3200];
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[ sample_num frame_num] = size(f);

set = zeros(frame_num, 12);

for k = 1 : frame_num
    window_frame = hanning_window(f(:,  k));
    fft_res = rfft(window_frame,fs,1);
    fft_res = abs(fft_res);
    
    for m = 1 : 12
        st = 0.0;
        for spectral = filter_bank(m) : filter_bank(m + 1) - 1
           st = st + fft_res(spectral);
        end
        set(k, m) = 10 * log10(st)/(filter_bank(m + 1) - filter_bank(m));
    end
    
end

end

function [ pitch_res ] = pitch( x, fs, frame_long, inc, alg_prog )
%short_energy_res, frame_long/2 * frame_num
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
p = mirpitch(f);
pitch_res = mirgetdata(p);
% pitch_res(isnan(pitch_res)) = 3000; 
% pitch_res = sort(pitch_res, 1);
[sample_num frame_num] = size(pitch_res);
pitch_res = pitch_res(1, :);

% for frame_index = 1 : frame_num
%     if(isnan(pitch_res(1, frame_index)))
%         pitch_res(1, frame_index) = 0;
%     end
% end
pitch_res = pitch_res';

end
% function [ pitch_res ] = pitch( x, fs, frame_long, inc, alg_prog )
% %short_energy_res, frame_long/2 * frame_num
% [fmap,pt2] = ftrack(x,fs);
% pitch_res = pt2';
% end

function [ spectrl_energy_res ] = spectrl_energy( x, fs, frame_long, inc, alg_prog )
%因为各篇论文中使用到的频率范围都不一致，所以本算法使用ERB带通能量
%mel_bank_filter_res frame_num * 15
a = miraudio(x, fs);
filter_res = mirfilterbank(a, 'Gammatone', 'NbChannels', 15);
f = mirframe(filter_res, frame_long, 'sp', inc,  'sp');
rms = mirrms(f);
rms = mirgetdata(rms);
spectrl_energy_res = [];
[sample_num frame_num filter_num] = size(rms);
temp = zeros(1, filter_num);
for frame_index = 1 : frame_num
    temp(1, :)  = rms(1, frame_index, :);
    spectrl_energy_res = [spectrl_energy_res;temp(1, :)];
end

end


function [ jitter_res ] = jitter( x, fs, frame_long, inc, alg_prog )
%jitter_res, frame_num * 1
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[sample_num frame_num] = size(f);
jitter_res = [];
for frame_index = 1 : frame_num
    [peaks peaks_pos] = find_peak_between_zeros(f(:, frame_index));
    [peaks_num dim] = size(peaks_pos);
    if(peaks_num > sample_num/50)
        T = zeros(peaks_num -1, 1);
        for peaks_index = 2 : peaks_num
            T(peaks_index - 1, 1) = peaks_pos(peaks_index, 1) - peaks_pos(peaks_index - 1, 1);
        end
        
        xxx = 0;
        yyy = 0;
        for T_index = 2 : peaks_num - 2
            xxx = xxx + 2 * T(T_index, 1) - T(T_index - 1, 1) - T(T_index + 1, 1);
            yyy = yyy + T(T_index, 1);
        end
        jitter_res = [jitter_res; xxx/yyy];
    else
        jitter_res = [jitter_res; 0];
    end
end

end

function [ shimmer_res ] = shimmer( x, fs, frame_long, inc, alg_prog )
%jitter_res, frame_num * 1,与JITTER唯一不一样的是，JITTER使用峰的间隔，而SHIMMER使用峰值
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[sample_num frame_num] = size(f);
shimmer_res = [];
for frame_index = 1 : frame_num
    [peaks peaks_pos] = find_peak_between_zeros(f(:, frame_index));
    [peaks_num dim] = size(peaks);
    if(peaks_num > sample_num/50)
        T = zeros(peaks_num -1, 1);
        for peaks_index = 2 : peaks_num
            T(peaks_index - 1, 1) = peaks(peaks_index, 1) - peaks(peaks_index - 1, 1);
        end
        
        xxx = 0;
        yyy = 0;
        for T_index = 2 : peaks_num - 2
            xxx = xxx + 2 * T(T_index, 1) - T(T_index - 1, 1) - T(T_index + 1, 1);
            yyy = yyy + T(T_index, 1);
        end
        shimmer_res = [shimmer_res; xxx/yyy];
    else
        shimmer_res = [shimmer_res; 0];
    end
end

end

function [ spectrum_res ] = mspectrum( x, fs, frame_long, inc, alg_prog )
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[sample_num frame_num] = size(f);
spectrum_res = [];
for frame_index = 1 : frame_num
    window_res = hanning_window(f(:, frame_index));
    c=cceps(window_res);
    c=fftshift(c);
    spectrum_res = [spectrum_res;c'];
end
end

function [ rspectrum_res ] = rspectrum( x, fs, frame_long, inc, alg_prog )
%x = awgn(x, 35);
% fft_x = fft(x);
% fft_x(length(fft_x)/8 : length(fft_x) * 7/8) = 0;
% ifft_x = ifft(fft_x);
% x = real(ifft_x);
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
[sample_num frame_num] = size(f);
rspectrum_res = [];
for frame_index = 1 : frame_num
     window_res = hanning_window(f(:, frame_index));
     d = rceps(window_res);
     d=fftshift(d);
%     d = mel_one_dim(window_res', fs);
    rspectrum_res = [rspectrum_res;d'];
end
 %rspectrum_res = rspectrum_preprocess(rspectrum_res);
end

function [ rspectrum_res ] = sun1_rspectrum( x, fs, frame_long, inc, alg_prog )
x = awgn(x, 35);
win=0.54+0.46*cos((1-frame_long:2:frame_long)*pi/frame_long);  % Hamming window
win=win/sqrt(sum(win.^2));      
[sf,t]=enframe(x,win,inc);
[frame_num feature_dim] = size(sf);
rspectrum_res = [];
for frame_index = 1 : frame_num
    d = rceps(sf(frame_index, :));
    rspectrum_res = [rspectrum_res;d(1 : floor((length(d) + 1)/2))];
end
 rspectrum_res = rspectrum_preprocess(rspectrum_res);
end

function [ short_fft_res ] = sun_rspectrum( x, fs, frame_long, inc, alg_prog )
lowest_freq = alg_prog;
[T,F,B] = spgrambw(x(:, 1),fs,'pjcwm', 60);
short_fft_res =log(B);
for index = 1 : length(F)
    if F(index) < lowest_freq 
        short_fft_res(:, index) = min(min(short_fft_res));
    end
end
end

function [ preprocess_res ] = rspectrum_preprocess(feature )
[frame_num feature_dim] = size(feature);
mean_features = mean(mean(feature));
feature(feature < mean_features) = mean_features;   
feature(:, 1 :floor(feature_dim/15)) = mean_features;
preprocess_res = feature;
end

function [ filter_res] = max_filter(features, window_size)
[rows cols] = size(features);
filter_res = zeros(rows, cols);

for row = floor(window_size/2) + 1 : rows - floor(window_size/2) - 1
    for col = floor(window_size/2) + 1: cols - floor(window_size/2) - 1
        winrow_start = row -  floor(window_size/2);
        winrow_end = row -  floor(window_size/2);
        wincol_start = col -  floor(window_size/2);
        wincol_end = col -  floor(window_size/2);    
        filter_res(row, col) = max(max(features(winrow_start : winrow_end, wincol_start :wincol_end )));
    end
end
end

function [ msf_res ] = msf(x, fs, frame_long, inc, alg_prog )
X = miraudio(x(:), fs);

filter_res = mirfilterbank(X, 'Bark');
f = mirframe(filter_res, frame_long, 'sp', inc,  'sp');
filter_res = mirgetdata(f );
[sample_num frame_num filter_bank_num] = size(filter_res);
H = zeros(5, filter_bank_num, frame_num);
for filter_bank_index = 1 : filter_bank_num
    for frame_index = 1 : 20:frame_num
        hin = hilbert_envelope(filter_res(:, frame_index, filter_bank_index));
        %t_x = miraudio(hin, fs);
        %t_x_filter_res = mirfilterbank(t_x, 'Gammatone', 'Manual', [-Inf 4 8 16 32 64 Inf]);
        %t_x_filter_res = mirgetdata(t_x_filter_res );
        fft_hin_res = fft(hin);
        fft_hin_res = abs(fft_hin_res(1 : floor(length(fft_hin_res)/2)));
    str_num = sprintf('%d_%d', filter_bank_index, frame_index);
    save_path = strcat('E:\datas_res\msf\', str_num, 'correct.jpg'); 
    h_fig = figure('Visible', 'off');
    hold on
    plot(fft_hin_res');  
    hold off
    saveas(h_fig, save_path);
    close(h_fig); %关闭figure，清空内存
%         for index = 1 : 5
%             sum_filter = sum(t_x_filter_res(:, :, index));
%             H(index,filter_bank_index, frame_index) = sum_filter;
%         end
    end
end
E = zeros(filter_bank_num, 5);
for filter_bank_index = 1 : filter_bank_num
    for index = 1 : 5
        E(filter_bank_index, index) = sum(H(index, filter_bank_index, :));
    end
end
Q1 = sum(E, 1)/filter_bank_num;
Q2 = zeros(filter_bank_num, 1);
for index = 1 : 5;
    Q2j = 1;
    for filter_bank_index = 1 : filter_bank_num
        Q2j = Q2j * E(filter_bank_index, index);
    end
    Q2(index) = Q2j^(1/filter_bank_num)/Q1(index);
end
Q4 = zeros(filter_bank_num, 1);
for index = 1 : 5;
    Q4lu = 0;
    Q4ld = 0;
    for filter_bank_index = 1 : filter_bank_num
        Q4lu = Q4lu + filter_bank_index * E(filter_bank_index, index);
        Q4ld = Q4ld + E(filter_bank_index, index);
    end
    Q4(index) = Q4lu/Q4ld;
end

msf_res = [Q1;Q2;Q4];
end


% function [ rspectrum_res ] = wave_trans( x, fs, frame_long, inc, alg_prog )
% a = miraudio(x, fs);
% f = mirframe(a, frame_long, 'sp', inc,  'sp');
% f = mirgetdata(f);
% [sample_num frame_num] = size(f);
% rspectrum_res = [];
% T = cell(frame_num, 1);
% fc = [];
% ac = []
% mc = [];
% for frame_index = 1 : frame_num
%     window_res = hanning_window(f(:, frame_index));
%     tree = wpdec(window_res,3,'db2');
%     Ti = [];
%     for coef_index = 0 : 4
%         X = wpcoef(tree,[3,coef_index]);
%         Tin = TEO(X');
%         Ti = [Ti;Tin];
%     end
%     T{frame_index, 1} = Ti;
%     mcn = mean(frame_T);
%     mc = [mc;mcn];
% end
% m = [];
% m = [m;1];
% for frame_index = 1 : frame_num
%    rn = T{frame_num - frame_index + 1, 1};
%    Tn = T{frame_index, 1};
%    fcn = 0;
%    fcd = 0;
%    acn = 0;
%    for coef_index = 1 : 5
%         rni = r(6 - i, :);
%         fcn = fcn + Tn(coef_index, :) * coef_index + rn(coef_index, :) * coef_index * m(frame_index);
%         fcd = fcd + Tn(coef_index, :) + rn(coef_index, :) * m(frame_index);
%         acn = fcn;
%    end 
%    fc = fc./ fcd;
%    ac = ac./55;
%    
%    if 
% end
% 
% end

function [ window_res ] = hanning_window( x )
%HANNING_WINDOW Summary of this function goes here
%   Detailed explanation goes here
[rows cols] = size(x);
h = hamming(rows);
window_res = zeros(rows, cols);
for n = 1 : cols
    for m = 1 : rows
        window_res(m, n) = h(m) * x(m, n);
    end
end
end

function [ Hkin] = hilbert_envelope( x )
%HILBERT_ENVELOPE Summary of this function goes here
%   Detailed explanation goes here
hil_res = hilbert(x);
Hkin = sqrt(imag(hil_res).^2 + real(hil_res).^2);
end

function [ peaks peaks_pos up_zerocross down_zerocross] = find_peak_between_zeros( x )
%FIND_PEAK_BETWEEN_ZEROS Summary of this function goes here
%   Detailed explanation goes here
[sample_num col] = size(x);
pre_sample_index = 1;
bk_sample_index = 2;
up_zerocross = [];
down_zerocross = [];
peaks = [];
peaks_pos = [];
max_peak = 0;
max_peak_pos = 0;
for sample_index = 2 : sample_num
    if x(sample_index, 1) == 0
        continue;
    end
    bk_sample_index = sample_index;
    if x(pre_sample_index, 1) <= 0 && x(bk_sample_index, 1) > 0
        up_zerocross = [up_zerocross; bk_sample_index];
        if max_peak > 0 %每次找到一个新的上升沿的零点时将这个上升沿与上一个上升沿之间的极大值增加到峰值数组中。此处的判断是为了跳过第一个峰值
            peaks = [peaks; max_peak];
            peaks_pos = [peaks_pos; max_peak_pos];
        end
        max_peak = x(bk_sample_index, 1);
        max_peak_pos = bk_sample_index;
    elseif x(pre_sample_index, 1) >= 0 && x(bk_sample_index, 1) < 0
        down_zerocross = [down_zerocross; bk_sample_index];
    end
    pre_sample_index = bk_sample_index;
    if max_peak > 0 && max_peak < x(bk_sample_index, 1) %当找到第一个上升沿的零点后，得到两上升沿之间的峰值
        max_peak = x(bk_sample_index, 1);
        max_peak_pos = bk_sample_index;
    end
end

end

function [ peaks  up_zerocross down_zerocross] = find_max_teo_between_zeros( x )
%FIND_PEAK_BETWEEN_ZEROS Summary of this function goes here
%   Detailed explanation goes here
[sample_num col] = size(x);
pre_sample_index = 1;
bk_sample_index = 2;
up_zerocross = [];
down_zerocross = [];
peaks = [];
max_peak = 0;
teo_res = TEO(x);
for sample_index = 2 : sample_num
    if x(sample_index, 1) == 0
        continue;
    end
    bk_sample_index = sample_index;
    if x(pre_sample_index, 1) <= 0 && x(bk_sample_index, 1) > 0
        up_zerocross = [up_zerocross; bk_sample_index];
        if max_peak > 0 %每次找到一个新的上升沿的零点时将这个上升沿与上一个上升沿之间的极大值增加到峰值数组中。此处的判断是为了跳过第一个峰值
            peaks = [peaks; max_peak];
        end
        max_peak = teo_res(bk_sample_index, 1);
    elseif x(pre_sample_index, 1) >= 0 && x(bk_sample_index, 1) < 0
        down_zerocross = [down_zerocross; bk_sample_index];
    end
    pre_sample_index = bk_sample_index;
    if max_peak > 0 && max_peak < teo_res(bk_sample_index, 1) %当找到第一个上升沿的零点后，得到两上升沿之间的峰值
        max_peak = teo_res(bk_sample_index, 1);
    end
end
end
