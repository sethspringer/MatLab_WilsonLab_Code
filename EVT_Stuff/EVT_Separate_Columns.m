function EVT_Separate_Columns

%Outputs two vector files, one for Tmu and the other for TriNo
%   Help save time in copying behavioral data for accuracy & RT
%   Seth Springer

[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

%create header for evt files%
triggers_header = 'TriNo';   
time_header = 'Tmu'; 
                                     
%set number of iterations based on whether there is one or multiple files%
cd(path);  
if isa(files,'char')                   %Checking to see if only one file is input (resulting in 'char' data type).                                                     
    iter = 1;
else
    iter = length(files);              %Setting # of iterations based on the # of selected files.
end

for i = 1:iter                                                                              
    if isa(files,'char')              %Checking to see if only one file is input (resulting in 'char' data type).                                                      
        data = readBESAevt(files);
    elseif length(files) > 1;
        data = readBESAevt(files{i});   %Setting variable "data" to be the ith file imported.
    else
        error('No compatible files selected! Boom roasted!')
    end
    
    %Creating separate vectors for the time and trigger columns.
    triggers = data(:,3);
    time = data(:,1);
    
    
    



%%SAVING FILES%%
%Creating a filename to which header and data will be saved.

        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_triggers.evt');   %Taking "file" and replacing ".evt" with "_recoded.evt".
        else
            filename = strcat(files{i}(1,1:end-4),'_triggers.evt');   %Taking the ith "file" and replacing ".evt" with "_recoded.evt".                       
        end
        
%Creating a blank file in the open directory with the "filename".
        fid = fopen(filename,'wt');
%Writing the .evt header to the blank file with name "filename".        
        fprintf(fid,'%s\n',triggers_header);        
        fclose(fid);
%Writing manipulated data, "evt_info", to the file with "filename" that already has a header.        
        dlmwrite(filename,triggers,'delimiter','\t','-append','precision','%.0f');
        
        
        
        
%Creating a filename to which timestamp header and data will be saved

if isa(files,'char')
      filename2 = strcat(files(1,1:end-4),'_Tmu.evt');   %Taking "file" and replacing ".evt" with "_recoded.evt".
else
      filename2 = strcat(files{i}(1,1:end-4),'_Tmu.evt');   %Taking the ith "file" and replacing ".evt" with "_recoded.evt".                       
end
        
%Creating a blank file in the open directory with the "filename2".
        fid2 = fopen(filename2,'wt');
%Writing the .evt header to the blank file with name "filename2".        
        fprintf(fid2,'%s\n',time_header);        
        fclose(fid2);
%Writing manipulated data, "time_when_trigger_is_12291_millisec", to the file with "filename2" that already has a header.        
        dlmwrite(filename2,time,'delimiter','\t','-append','precision','%.0f');    
end

disp('DONE!')
end

