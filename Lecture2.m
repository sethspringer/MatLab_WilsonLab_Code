function Lecture2


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

%create header for evt files%
evt_header = 'Tmu         	Code	TriNo';   
evt_header_timestamp = 'Time (ms)'; 
                                     
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
    
    
    
%%%%Actual manipulations to data%%%%
    %Timestamping all instances of 12500s.
    time_when_trigger_is_12500s_microsec = time(triggers > 12499 & triggers < 12600);
    time_when_trigger_is_12500s_millisec = time_when_trigger_is_12500s_microsec./1000; %Conversion
    
    %Timestamping all instances of 12800s.
    time_when_trigger_is_12800s_microsec = time(triggers > 12799 & triggers < 12900);
    time_when_trigger_is_12800s_millisec = time_when_trigger_is_12800s_microsec./1000; %Conversion
    
    %Converting 12289, 12290, 12291, and 12292 to 12290.
    triggers(triggers < 12293 & triggers > 12288) = 12290;
    
    
evt_info = [time,data(1:size(time,1),2),triggers];    %Concatenating separate .evt vectors.



%%SAVING FILES%%
%Creating a filename to which header and data will be saved.

        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_recoded_for_lecture2.evt');   %Taking "file" and replacing ".evt" with "_recoded.evt".
        else
            filename = strcat(files{i}(1,1:end-4),'_recoded_for_lecture2.evt');   %Taking the ith "file" and replacing ".evt" with "_recoded.evt".                       
        end
        
%Creating a blank file in the open directory with the "filename".
        fid = fopen(filename,'wt');
%Writing the .evt header to the blank file with name "filename".        
        fprintf(fid,'%s\n',evt_header);        
        fclose(fid);
%Writing manipulated data, "evt_info", to the file with "filename" that already has a header.        
        dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
        
        
        
        
%Creating a filename to which 12500s timestamp header and data will be saved

if isa(files,'char')
      filename2 = strcat(files(1,1:end-4),'_timestamp_12500s.evt');   %Taking "file" and replacing ".evt" with "_recoded.evt".
else
      filename2 = strcat(files{i}(1,1:end-4),'_timestamp_12500s.evt');   %Taking the ith "file" and replacing ".evt" with "_recoded.evt".                       
end
        
%Creating a blank file in the open directory with the "filename2".
        fid2 = fopen(filename2,'wt');
%Writing the .evt header to the blank file with name "filename2".        
        fprintf(fid2,'%s\n',evt_header_timestamp);        
        fclose(fid2);
%Writing manipulated data, "time_when_trigger_is_12500s_millisec", to the file with "filename2" that already has a header.        
        dlmwrite(filename2,time_when_trigger_is_12500s_millisec,'delimiter','\t','-append','precision','%.0f');    

        
%Creating a filename to which 12800s timestamp header and data will be saved

if isa(files,'char')
      filename3 = strcat(files(1,1:end-4),'_timestamp_12800s.evt');   %Taking "file" and replacing ".evt" with "_recoded.evt".
else
      filename3 = strcat(files{i}(1,1:end-4),'_timestamp_12800s.evt');   %Taking the ith "file" and replacing ".evt" with "_recoded.evt".                       
end
        
%Creating a blank file in the open directory with the "filename2".
        fid3 = fopen(filename3,'wt');
%Writing the .evt header to the blank file with name "filename2".        
        fprintf(fid3,'%s\n',evt_header_timestamp);        
        fclose(fid3);
%Writing manipulated data, "time_when_trigger_is_12800s_millisec", to the file with "filename2" that already has a header.        
        dlmwrite(filename3,time_when_trigger_is_12800s_millisec,'delimiter','\t','-append','precision','%.0f');    
        
        
        
end

disp('DONE!')



end

