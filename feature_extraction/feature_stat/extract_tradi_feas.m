function [ features ] = extract_tradi_feas(  x, fs, frame_long, inc, alg_type, alg_prog )
%EXTRACT_TRADI_FEAS Summary of this function goes here
%   Detailed explanation goes here
%alg_type, �㷨����, Ŀǰ�У�'short_energy', 'short_zerocross', 'short_correct',
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
        %short_fft_res, frame_num * frame_long������С��2��N�η�+1
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
        %frame_num * 12, 12ά�˲�����
        features = LFPC(x,  fs, frame_long, inc, alg_prog);
    case 'sun_mfcc'
        %frame_num * alg_prog
        features = sun_mfcc(x, fs, frame_long, inc, alg_prog);
    case 'sun_lfpc'
        %frame_num * alg_prog
        features = sun_lfpc(x, fs, frame_long, inc, alg_prog);        
    case 'autocorr_contours'
         %frame * 1
        %����������VAD,Ҳ��Ŀǰ�ռ�������õ�VAD�㷨������0.5Ϊ�����Ƭ��,���Ǹ÷���ֻ�ܶ�����������Ч
        %x = awgn(x, 20);
        features = autocorr_contours(x, fs, frame_long, inc, alg_prog);
    case 'harmfreq_MOLRT_VAD'
        %���VADЧ���ܲ�
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
        %��������������ʱû������֤
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
%short_fft_res, frame_long������С��2��N�η�+1 * frame_num
 [T,F,B] = spgrambw(x(:, 1),fs,'pjcwm', 60);
 short_fft_res =log10(B);
end
%lpcc�������ϵ�������
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
ar = ar(1, 2 : p);%��Ϊ����LPC�Ŀ���ar��һ��λ������Ϊ�������1������ȥ��
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
ar = ar(1, 2 : p);%��Ϊ����LPC�Ŀ���ar��һ��λ������Ϊ�������1������ȥ��
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
%��Ϊ��Ҫ�����ʣ�����X��matlab�Դ���wavread��ȡ��
%ZCPA������
%��һ������֡
%�ڶ������˲�
%����������ÿ���˲��������ÿ��֡�����ʱ������
%���Ĳ���ͨ�������ʵõ����
%���岽������ֱ��ͼ
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
%��һ������ÿһ֡��amplitude value
for t = 1 : size(f, 2)
    ft = f(:, t);
    A = [A;sum(ft.^2)]; 
end
A(vad < 0.5) = NaN;
%�ڶ�������ƽ��amplitude value
mean_A = mean(A(~isnan(A)));
%������:��shimmer
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
        if max_peak > 0 %ÿ���ҵ�һ���µ������ص����ʱ���������������һ��������֮��ļ���ֵ���ӵ���ֵ�����С��˴����ж���Ϊ��������һ����ֵ
            peaks = [peaks; max_peak];
            peaks_pos = [peaks_pos; max_peak_pos];
        end
        max_peak = x(bk_sample_index, 1);
        max_peak_pos = bk_sample_index;
    elseif x(pre_sample_index, 1) >= 0 && x(bk_sample_index, 1) < 0
        down_zerocross = [down_zerocross; bk_sample_index];
    end
    pre_sample_index = bk_sample_index;
    if max_peak > 0 && max_peak < x(bk_sample_index, 1) %���ҵ���һ�������ص����󣬵õ���������֮��ķ�ֵ
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
        if max_peak > 0 %ÿ���ҵ�һ���µ������ص����ʱ���������������һ��������֮��ļ���ֵ���ӵ���ֵ�����С��˴����ж���Ϊ��������һ����ֵ
            peaks = [peaks; max_peak];
        end
        max_peak = teo_res(bk_sample_index, 1);
    elseif x(pre_sample_index, 1) >= 0 && x(bk_sample_index, 1) < 0
        down_zerocross = [down_zerocross; bk_sample_index];
    end
    pre_sample_index = bk_sample_index;
    if max_peak > 0 && max_peak < teo_res(bk_sample_index, 1) %���ҵ���һ�������ص����󣬵õ���������֮��ķ�ֵ
        max_peak = teo_res(bk_sample_index, 1);
    end
end
end
