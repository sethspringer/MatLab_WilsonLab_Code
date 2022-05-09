function seth_preproc_moveIntoSubjectFolders_pullBy_Name_Data

%PURPOSE:           Move files into folders with the same starting name (MAKES FOLDERS IF THEY DON'T EXIST).
%
%REQUIRED INPUTS:   FIF files.
%
%		    
%NOTES:            
%		                
%AUTHOR:            Seth D. Springer, DICoN Lab, Boys Town National Research Hospital
%VERSION HISTORY:   2/3/2022  v1: First working version of program
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[files,path,~] = uigetfile('*.fif','Please Select FIFs','Multiselect','on');        %select evt files%

cd(path);


n = length(files);

for i = 1:n
    
    guid = files{i}(1:9);
    guid_wild = strcat(guid,'*');
    
    scan_date = files{i}(16:23);
    
    scan_date_num = str2num(scan_date);
    
    
    
    if scan_date_num > 20211101
       movefile(guid_wild,'FIFs');
    end
    
    
end %end of subject for loop 





end %end of function
