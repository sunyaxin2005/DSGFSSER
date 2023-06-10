function [ train_inxs ] = get_train_inxs( cross_type, emotion_label, people_label )
%GET_TRAIN_INXS Summary of this function goes here
%   Detailed explanation goes here
t=emotion_label-min(emotion_label)+1;
groups =t'; 
switch cross_type
    case 'ten_cross'
        rand('state',sum(100*clock)*rand(1));
        indices = crossvalind('Kfold',groups,10);
    case 'five_cross'
        rand('state',sum(100*clock)*rand(1));
        indices = crossvalind('Kfold',groups,5);
    case 'LOSO'
        indices = people_label;
    case 'LOGSO'
        indices = people_label;
        train_inxs = cell(5, 1);
          male = [1, 4, 5, 6, 9];
          fmale = [2, 3, 7, 8, 10];
          for j = 1 : 5
            test = (indices == male(5 - j + 1) | indices == fmale(j)); train = ~test;
            train_inxs{j} = train;
          end
end
uinx = unique(indices);
train_inxs = cell(length(uinx), 1);
for inx = 1 : length(train_inxs)
    test = (indices == uinx(inx));
    train = ~test;
    train_inxs{inx} = train;
end

end

