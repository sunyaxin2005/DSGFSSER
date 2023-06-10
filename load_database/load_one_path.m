function [ files_name ] = load_one_path(  wave_files_names, wave_files_path   )
%LOAD_ONE_PATH Summary of this function goes here
%   Detailed explanation goes here
file_num = length(wave_files_names)
files_name = cell(file_num, 1);
% for ci = 1 : file_num
%     files_name{ci} = wave_files_names(ci).name;
% end
% files_name = sort_string(files_name);
for ci = 1 : file_num
    files_name{ci} = strcat(wave_files_path,  wave_files_names(ci).name);
end

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