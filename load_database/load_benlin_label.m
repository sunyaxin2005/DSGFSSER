function [ files_path files_people_label files_emotion_label files_sentence_label words_mark ] = load_benlin_label( wave_files_names, wave_files_path  )
%LOAD_BENLIN_LABEL Summary of this function goes here
%   Detailed explanation goes here
load_type = 'load_all';
file_num = length(wave_files_names)
files_name = cell(file_num, 1);
for ci = 1 : file_num
    files_name{ci} = wave_files_names(ci).name;
end
files_name = sort_string(files_name);
switch load_type 
    case 'load_all'
        [ files_path files_people_label files_emotion_label files_sentence_label words_mark] = load_all( files_name, wave_files_path );
    case 'load_same_sentece'
    case 'load_same_emotion'
    case 'load_same_people'
end
end

function [ files_path files_people_label files_emotion_label files_sentence_label words_mark] = load_all( files_name, wave_files_path )
load_file_num = length(files_name);
files_path = cell(load_file_num, 1);
files_people_label = cell(load_file_num, 1);
files_emotion_label = cell(load_file_num, 1);
files_sentence_label = cell(load_file_num, 1);
index = 1;
words_mark = cell(load_file_num, 1);
for ci = 1:1:load_file_num
    file_path = strcat(wave_files_path, files_name{ci});
    files_path{index} = file_path;
    file_name = files_name{ci};
    file_name_length = length(file_name);
    files_people_label{index} = file_name(1:2);
    files_emotion_label{index} = file_name(6);
    files_sentence_label{index} = file_name(3:5);
%     [word_mark words_content] = load_word_mark( file_path );    
%      words_mark(index) = {word_mark}; 
    index = index + 1;
end
[ files_people_label files_emotion_label files_sentence_label]  = str_label2_num_label(files_people_label, files_emotion_label, files_sentence_label);
end

function [ files_people_label_num files_emotion_label_num files_sentence_label_num] = str_label2_num_label(files_people_label_str, files_emotion_label_str, files_sentence_label_str)
files_num = length(files_people_label_str);
files_people_label_num = zeros(files_num, 1);
files_emotion_label_num = zeros(files_num, 1);
files_sentence_label_num = zeros(files_num, 1);
people_have_strlabel = cell(50, 1);
emotion_have_strlabel = cell(50, 1);
sentence_have_strlabel = cell(50, 1);
people_num = 0;
emotion_num = 0;
sentence_num = 0;
for files_index = 1 : files_num
    [is_in label_index] = is_in_label(files_people_label_str{files_index}, people_have_strlabel, people_num);
    if is_in == 1
        files_people_label_num(files_index, 1) = label_index;
    else
        people_num = people_num + 1;
        files_people_label_num(files_index, 1) = people_num;
        people_have_strlabel{people_num} = files_people_label_str{files_index};
    end
    
    [is_in label_index] = is_in_label(files_emotion_label_str{files_index}, emotion_have_strlabel, emotion_num);
    if is_in == 1
        files_emotion_label_num(files_index, 1) = label_index;
    else
        emotion_num = emotion_num + 1;
        files_emotion_label_num(files_index, 1) = emotion_num;
        emotion_have_strlabel{emotion_num} = files_emotion_label_str{files_index};
    end
    
    [is_in label_index] = is_in_label(files_sentence_label_str{files_index}, sentence_have_strlabel, sentence_num);
    if is_in == 1
        files_sentence_label_num(files_index, 1) = label_index;
    else
        sentence_num = sentence_num + 1;
        files_sentence_label_num(files_index, 1) = sentence_num;
        sentence_have_strlabel{sentence_num} = files_sentence_label_str{files_index};
    end
end
end

function [is_in label_index]= is_in_label(str, label, num)
for index = 1 : num
    if strcmp(str, label{index}) == 1
       is_in = 1;
       label_index = index;
       return ;
    end
end
is_in = 0;
label_index = -1;
end

function [files_name] =  sort_string(files_name)
file_num = length(files_name);
for indexi = 1 :  file_num - 1
    for indexj = indexi + 1 : file_num
        if strm( files_name{indexi}, files_name{indexj}) > 0
           temp = files_name{indexi};
           files_name{indexi} = files_name{indexj};
           files_name{indexj} = temp;
        end
    end
end
end

function [p]=strm(str1,str2)
k=min(length(str1),length(str2));
for n=1:k   %比较前k位    
    if(str1(n)>str2(n))        
        p=1;break;    
    elseif(str1(n)==str2(n))       
        p=0;    
    else
        p=-1;break;    
    end
end
if(p==0)    
    if(length(str1)>length(str2)) %前k位相等，但str1更长       
        p=1;    
    elseif(length(str1)==length(str2))        
        p=0;   
    else
        p=-1;
    end
end
end