function [ features ] = extract_tradi_feas(  x, fs, frame_long, inc, alg_type, alg_prog )
%EXTRACT_TRADI_FEAS Summary of this function goes here
%   Detailed explanation goes here
%alg_type, 算法类型, 目前有：'short_energy', 'short_zerocross', 'short_correct',
%'short_fft', 'LPCC', 'MFCC', 'ZCPA', 'LFPC', 'sun_mfcc', 'sun_lfpc',
%'autocorr_contours', 'harmfreq_MOLRT_VAD', 'best_pitch', 'plp', 'RASTA-PLP', 'fromant',
%'jitter', 'shimmer'
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
    case 'LPCC'
        %frame_num * p
        features = LPC(x, fs, frame_long, inc, 'LPCC');
    case 'MFCC'
        %frame_num * 13
        features = mfcc(x, fs, frame_long, inc, alg_prog);
    case 'ZCPA'
        %frame_num * 16
        features = ZCPA( x,  fs, frame_long, inc, 10); 
    case 'LFPC'
        %frame_num * 12, 12维滤波器数
        features = LFPC(x,  fs, frame_long, inc, alg_prog);
    case 'sun_mfcc'
        %frame_num * alg_prog
        features = sun_mfcc(x, fs, frame_long, inc, alg_prog);
    case 'sun_lfpc'
        %frame_num * alg_prog
        features = sun_lfpc(x, fs, frame_long, inc, alg_prog);        
    case 'autocorr_contours'
         %frame * 1
        %可以拿来做VAD,也是目前收集到的最好的VAD算法，大于0.5为活动语音片段,但是该方法只能对浊音部分有效
        %x = awgn(x, 20);
        features = autocorr_contours(x, fs, frame_long, inc, alg_prog);
    case 'harmfreq_MOLRT_VAD'
        %这个VAD效果很差
        [ features ] = VAD_main( x, fs );
    case 'best_pitch'
         %frame * 1
        [features, score] = fast_mbsc_fixedWinlen_tracking(x,fs);
    case 'plp'
        % Calculate 12th order PLP features without RASTA
        [cep2, features] = rastaplp(x, fs, 1, 12);
        features = 10*log10(features);
        features = features';
        %imagesc(10*log10(features));
    case 'RASTA-PLP'
        %frame* 
        % Calculate basic RASTA-PLP cepstra and spectra
        [cep1, features] = rastaplp(x, fs, 1, 12); 
        features = 10*log10(features);
        features = features';
    case 'fromant'
        %frame * 3
        [features,pt2] = ftrack(x,fs);
        %以下两个特征暂时没经过验证
    case 'jitter'
        %frame * 1
        [ features ] = jitter(x, fs);
    case 'shimmer'
        [ features ] = shimmer(x, fs, frame_long, inc, alg_prog);
    case 'Q'%the default w = 5
        features = f_Q(x, fs, frame_long, inc, alg_prog); 
    case 'TQ'%the default w = 5
        features = f_TQ(x, fs, frame_long, inc, alg_prog); 
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
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
correct_res = mirautocor(f);
short_correct_res = mirgetdata(correct_res);
short_correct_res = short_correct_res';
end

function [ mfcc_res ] = mfcc( x, fs, frame_long, inc, alg_prog )
%mfcc, frame_num * mfcc_coefs
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
mfcc_res = mirmfcc(f);
mfcc_res = mirgetdata(mfcc_res);
mfcc_res = mfcc_res';
mfcc_res = mfcc_res(:, 2:13);
end

function [ short_fft_res ] = short_fft( x, fs, frame_long, inc, alg_prog )
%short_fft_res, frame_long向下最小的2的N次方+1 * frame_num
 [T,F,B] = spgrambw(x(:, 1),fs,'pjcwm', 60);
 short_fft_res =log10(B);
end
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
    ar = lpcauto(f(:, frame_index), 14);
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

function [ b ] = sun_mfcc( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2: 13)];
end
b = mfcc_res;
end

function [ b ] = sun_lfpc( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
end

function [ b ] = f_Q( x, fs, frame_long, inc, alg_prog )
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = speech_HU( b );
b = b1;
end
function [ b ] = f_TQ( x, fs, frame_long, inc, alg_prog )
%x = awgn(x, 30);
b = mel_fft( x, fs, frame_long, inc, alg_prog);
[b1, b2, b3, b4, b5, b6, b7] = speech_HU( b );
b = b1;
mfcc_res = [];
for frame_index = 1 : size(b, 1)
    dct_res = dct(b(frame_index, :));
    mfcc_res = [mfcc_res; dct_res(2 : 13)];
end
b = mfcc_res;
end

function [ features ] = jitter(x, fs, frame_long, inc, alg_prog )
[pitch, score] = fast_mbsc_fixedWinlen_tracking(x,fs);
pitch(pitch == 0) = NaN;
mean_pitch = sum(pitch(~isnan(pitch)))/length(pitch(~isnan(pitch)));
features = zeros(length(pitch), 1);
for inx = 1 : length(pitch) - 1
    features(inx) = abs(pitch(inx) - pitch(inx + 1))/mean_pitch;
end
features(length(pitch)) = NaN;
end

function [ features ] = shimmer(x, fs, frame_long, inc, alg_prog )
a = miraudio(x, fs);
f = mirframe(a, frame_long, 'sp', inc,  'sp');
f = mirgetdata(f);
A = [];
vad = autocorr_contours(x, fs, frame_long, inc, alg_prog);
%第一步：求每一帧的amplitude value
for t = 1 : size(f, 2)
    ft = f(:, t);
    A = [A;sum(ft.^2)]; 
end
A(vad < 0.5) = NaN;
%第二步：求平均amplitude value
mean_A = mean(A(~isnan(A)));
%第三步:求shimmer
features = zeros(length(A), 1);
for inx = 1 : size(f, 2) - 1
    features(inx) = abs(A(inx) - A(inx + 1))/mean_A;
end
features(length(A)) = NaN;
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
