function [ output_args ] = draw_plot( features )
%DRAW_PLOT Summary of this function goes here
%   Detailed explanation goes here
[frame_num features_num] = size(features);
h_fig = figure('Visible', 'on');
x = sqrt(features_num);
x1 = floor(x);
if(x > x1)
    x = x1 + 1;
end
linspec{1, 1} = '--*r';
linspec{2, 1} = '-.*g';
linspec{3, 1} = '-.*b';
linspec{4, 1} = '-.*c';
linspec{5, 1} = '-.*m';
linspec{6, 1} = '-.*k';
linspec{7, 1} = ':>b';
linspec{8, 1} = ':>g';
linspec{9, 1} = ':>b';
linspec{10, 1} = ':>c';
linspec{11, 1} = ':>m';
linspec{12, 1} = ':>k';
linspec{13, 1} = ':or';
linspec{14, 1} = ':og';
linspec{15, 1} = ':ob';
linspec{16, 1} = ':oc';
linspec{17, 1} = ':om';
linspec{18, 1} = ':ok';
if features_num < 19
    hold on
    x = (10:10:450);
    for features_index = 1 : features_num
    %subplot(x, x, features_index);
        plot(x, features(:, features_index), linspec{features_index});
    end
    %legend('SSMCFS', 'SMCFS', 'MCFS', 'PCA', 'SFS', 'mRMR');
    %legend('SSMCFS', 'PCA', 'mRMR', 'LS', 'MCFS', 'DISR', 'SFS');
    legend('PROS+HU', 'PROS+LPCC', 'PROS+MFCC', 'PROS+PLP', 'PROS', 'PROS+ZCPA');
    xlabel('# of features');
    ylabel('precision')
    %saveas(h_fig, save_path);
    title('Speaker Independent And Speaker Normal(Berlin)')
    %axis([10 350 0.65 0.89])
        Aaa = legend('ESELF', 'SDA', 'LRGPSSDR', 'PCA', 'LDA', 'LRLFSDR');
    set(Aaa,'FontSize', 12);
    Aaa = xlabel('# of features');
    set(Aaa,'FontSize', 12);
    Aaa = ylabel('precision');
    set(Aaa,'FontSize', 12);
    %saveas(h_fig, save_path);
    Aaa = title('Speaker Independent And Speaker Normal(Berlin)');
    set(Aaa,'FontSize', 12);
    
    hold off
% else 
%     for features_index = 1 : features_num
%         subplot(x, x, features_index);
%         plot(features(:, features_index));
%     end
%     saveas(h_fig, save_path);
% end
%close(h_fig); %¹Ø±Õfigure£¬Çå¿ÕÄÚ´æ 
end

end

