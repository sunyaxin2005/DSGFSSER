function [ output_args ] = draw_plot_DR( features )
%DRAW_PLOT Summary of this function goes here
%   Detailed explanation goes here
[frame_num features_num] = size(features);
h_fig = figure('Visible', 'on');
x = sqrt(features_num);
x1 = floor(x);
if(x > x1)
    x = x1 + 1;
end
% linspec{6, 1} = '--*r';
% linspec{2, 1} = '-.*g';
% linspec{3, 1} = '-.*b';
% linspec{4, 1} = '-.*c';
% linspec{5, 1} = '-.*m';
% linspec{1, 1} = '-.*k';
% linspec{7, 1} = ':>b';
% linspec{8, 1} = ':>g';
% linspec{9, 1} = ':>b';
% linspec{10, 1} = ':>c';
% linspec{11, 1} = ':>m';
% linspec{12, 1} = ':>k';
% linspec{13, 1} = ':or';
% linspec{14, 1} = ':og';
% linspec{15, 1} = ':ob';
% linspec{16, 1} = ':oc';
% linspec{17, 1} = ':om';
% linspec{18, 1} = ':ok';
% linspec{1, 1} = '-ok';
% linspec{2, 1} = '--*k';
% linspec{3, 1} = '--+k';
% linspec{4, 1} = '--.k';
% linspec{5, 1} = '--xk';
% linspec{6, 1} = '-->k';
% linspec{7, 1} = '--<k';
% linspec{8, 1} = ':>k';
% linspec{9, 1} = ':>k';

linspec{1, 1} = '-*r';
linspec{2, 1} = '-*g';
linspec{3, 1} = '-*b';
linspec{4, 1} = '-*c';
linspec{5, 1} = '-*m';
linspec{6, 1} = '-*k';
if features_num < 19
    hold on
    %x = (10:10:450);
    for features_index = 1 : features_num
        temp_fea = features(:, features_index);
%         len = length(temp_fea(temp_fea > 0.1) );
%         temp_fea = temp_fea(1:len);
%         x = (10:20:20 * len);
max_dim = 140;
    x = 10:10:140;
    %subplot(x, x, features_index);
        plot(x,temp_fea, linspec{features_index});
    end
    
    set(gca,'FontSize', 12);
    %legend('SSMCFS', 'SMCFS', 'MCFS', 'PCA', 'SFS', 'mRMR');
    %legend('SSMCFS', 'PCA', 'mRMR', 'LS', 'MCFS', 'disr');
    %Aaa = legend('ESELF', 'SDA', 'LRGPSSDR', 'PCA', 'LDA', 'LRLFSDR');
    %Aaa = legend('SPFS', 'MRSSPFS6', 'MRSSPFS7', 'MRSSPFS8', 'MRSSPFS9');
    %Aaa = legend('AWSRC','WSRC','SRC','SVM','KNN');
    Aaa = legend('LRLFDSDR','PCA','NPE','LPP');
    %Aaa = legend('NN', 'SVM', 'SRC', 'SPPDSRC', 'WSRC', 'LSSRC');
    %Aaa = legend('MRSSPFS', 'MCFS', 'conred', 'MRMR', 'DISR','JMI', 'CMIFS');
    set(Aaa,'FontSize', 12);
    %Aaa = xlabel('# of features');
    Aaa = xlabel('Dimension');
    set(Aaa,'FontSize', 12);
    Aaa = ylabel('Accuracy');
    set(Aaa,'FontSize', 12);
    %saveas(h_fig, save_path);
    %Aaa = title('Speaker Independent And Speaker Normal(Berlin)');
    Aaa = title('Speaker Dependent(EmoDB)');
    set(Aaa,'FontSize', 12);
    axis([10 max_dim min(min(features))-0.01 max(max(features)) + 0.01])
    xlim([10, max_dim]);
    set(gca,'XTick', [20:20:max_dim]);
    set(gca,'YTick', [0.62:0.04:0.78]);
    ylim([0.62,0.78]);
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

