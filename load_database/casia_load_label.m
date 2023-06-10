function [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = casia_load_label( wave_files_names, wave_files_path )
%CASIA_LOAD_LABEL Summary of this function goes here
%   Detailed explanation goes here
load_type = 'load_all';
switch load_type 
    case 'load_all'
        [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_all( wave_files_names, wave_files_path );
    case 'load_same_sentece'
        [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_sentece( wave_files_names, wave_files_path);
    case 'load_same_emotion'
        [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_emotion( wave_files_names, wave_files_path);
    case 'load_same_people'
        [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_people( wave_files_names, wave_files_path);
end


end

function [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_all( wave_files_names, wave_files_path )
[rows, cols] = size(wave_files_names);
files_name = cell(1200, 1);
files_people_label = zeros(1200, 1);
files_emotion_label = zeros(1200, 1);
files_sentence_label = zeros(1200, 1);
index = 1;
words_mark = cell(1200, 1);
for ci = 1:1:1200
    file_path = strcat(wave_files_path, wave_files_names(ci).name);
    files_name{index} = file_path;
    files_people_label(index) = floor((ci - 1)/300) + 1;
    files_emotion_label(index) = floor(rem(ci - 1,300)/50) + 1;
    files_sentence_label(index) = rem(ci - 1, 50) + 1;
    [word_mark words_content] = load_word_mark( file_path );    
     words_mark(index) = {word_mark}; 
    index = index + 1;
end

end

function [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_sentece( wave_files_names, wave_files_path )
[rows, cols] = size(wave_files_names);
files_name = cell(24, 1);
files_people_label = zeros(24, 1);
files_emotion_label = zeros(24, 1);
files_sentence_label = zeros(24, 1);
index = 1;
words_mark = cell(24, 1);
for ci = 1:50:1200
    file_path = strcat(wave_files_path, wave_files_names(ci).name);
    files_name{index} = file_path;
    files_people_label(index) = floor(ci/300) + 1;
    files_emotion_label(index) = floor(rem(ci - 1,300)/50) + 1;
    files_sentence_label(index) = rem(ci - 1, 50) + 1;
    [word_mark words_content] = load_word_mark( file_path );    
     words_mark(index) = {word_mark}; 
    index = index + 1;
end
end

function [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_emotion( wave_files_names, wave_files_path )

[rows, cols] = size(wave_files_names);
files_name = cell(200, 1);
files_people_label = zeros(200, 1);
files_emotion_label = zeros(200, 1);
files_sentence_label = zeros(200, 1);
index = 1;
words_mark = cell(200, 1);
for people_begin = 0:300:1199
for ci = people_begin + 1:1: people_begin + 50
    file_path = strcat(wave_files_path, wave_files_names(ci).name);
    files_name{index} = file_path;
    files_people_label(index) = floor(ci/300) + 1;
    files_emotion_label(index) = floor(rem(ci - 1,300)/50) + 1;
    files_sentence_label(index) = rem(ci - 1, 50) + 1;
    [word_mark words_content] = load_word_mark( file_path );    
     words_mark(index) = {word_mark}; 
    index = index + 1;
end
end

end

function [ files_name files_people_label files_emotion_label files_sentence_label words_mark] = load_same_people( wave_files_names, wave_files_path )

[rows, cols] = size(wave_files_names);
files_name = cell(60, 1);
files_people_label = zeros(60, 1);
files_emotion_label = zeros(60, 1);
files_sentence_label = zeros(60, 1);
index = 1;
words_mark = cell(60, 1);
for ci = 1:5:300
    file_path = strcat(wave_files_path, wave_files_names(ci).name);
    files_name{index} = file_path;
    files_people_label(index) = floor(ci/300) + 1;
    files_emotion_label(index) = floor(rem(ci - 1,300)/50) + 1;
    files_sentence_label(index) = rem(ci - 1, 50) + 1;
    [word_mark words_content] = load_word_mark( file_path );    
     words_mark(index) = {word_mark}; 
    index = index + 1;
end

end

function [words_mark words_content] = load_word_mark( file_path )
ind_dot = find(file_path == '.', 1, 'last');
file_path = file_path(1 : ind_dot);
peak_path = strcat(file_path, 'peak');
[a1, a2, a3] =  textread(peak_path, '%d %d %s');

words_num = length(a1);
words_mark = [];
words_content = {};

for words_index = 1 : words_num
    word = a3(words_index, 1);
    if strcmp(word, 'SIL') == 0
        word_mark = [a1(words_index), a2(words_index)];
        words_mark = [words_mark;word_mark];
        words_content = {words_content; a3(words_index)};
    end
end

end
%
% [rows, cols] = size(wave_files_names);
% files_name = cell(200, 1);
% files_people_label = zeros(200, 1);
% files_emotion_label = zeros(200, 1);
% files_sentence_label = zeros(200, 1);
% index = 1;
% for people_begin = 0:300:1199
% for ci = people_begin + 1:1: people_begin + 50
%     file_path = strcat(wave_files_path, wave_files_names(ci).name);
%     files_name{index} = file_path;
%     files_people_label(index) = floor(ci/300) + 1;
%     files_emotion_label(index) = floor(rem(ci - 1,300)/50) + 1;
%     files_sentence_label(index) = rem(ci - 1, 50) + 1;
%     index = index + 1;
% end
% end