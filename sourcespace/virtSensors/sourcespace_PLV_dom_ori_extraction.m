function sourcespace_PLV_dom_ori_extraction


%PURPOSE: Process PLV files by finding dominant orientaion for each
%participant depending on the methof of choice (peak max or average in time
%of interest)
%
%REQUIRED INPUTS:   .TFCs files from PLV batch in BESA.
%
%		    
%
%NOTES:            
%		   
%
%
%                  
%AUTHOR:   Camilo A. Castelblanco and Seth D. Springer, DICoN Lab, Boys Town National Research Hospital         
%VERSION HISTORY:  4/26/2022   v1: First working version of program
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prompt = {'Start Frequency','End Frequency'};
title = 'A time series will be extracted, enter the start and end frequency';
dims = [1 90];
definput = {'0','0'};
answer = inputdlg(prompt,title,dims,definput);
[start_freq,end_freq] = deal(answer{:});
start_freq =str2num(start_freq);
end_freq =str2num(end_freq);

if start_freq >= end_freq
    error('enter proper frequencies')
end


%Load files for source 1 orientation 1
[files_R11,path_R11,~] = uigetfile('*.tfc','Please Select The TFCs for Source 1 Orientation 1 (R11)','Multiselect','on');        %select evt files%



%Load files for source 1 orientation 2
[files_R12,path_R12,~] = uigetfile('*.tfc','Please Select The TFCs for Source 1 Orientation 2 (R12)','Multiselect','on');        %select evt files%


n_R11 = length(files_R11);
n_R12 = length(files_R12);

if n_R11 ~= n_R12
    error('The number of files for each orientation needs to be equal... Better luck next time...')
end



for ii = 1:2
    
    if ii == 1
        cd(path_R11)
    elseif ii == 2
        cd(path_R12)
    end
    
    for i=1:n_R11
        
        
        %---Create Average Waveform For Each File---%
        
        if ii == 1
            Str=fileread(files_R11{i});
            fid=fopen(files_R11{i},'rt');
            files = files_R11;
            
        elseif ii == 2
            Str=fileread(files_R12{i});
            fid=fopen(files_R12{i},'rt');
            files = files_R11;
            
        end
        
        %-----------------------------------------------------------------%
        %                     GET RID OF THE HEADER                       %
        %-----------------------------------------------------------------%
        headerLine1 = fgetl(fid);
        %disp(headerLine1)
        headerLine2 = fgetl(fid);
        %disp(headerLine2)
        
        
        %---Find file parameters---%
        Str(strfind(Str,'='))=[];
        
        Key='NumberTimeSamples';
        Index=strfind(Str,Key);
        numTimeSamples(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='FreqIntervalInHz';
        Index=strfind(Str,Key);
        FreqInterval(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='NumberFrequencies';
        Index=strfind(Str,Key);
        numFrequencies(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='TimeStartInMS';
        Index=strfind(Str,Key);
        startTime(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='IntervalInMS';
        Index=strfind(Str,Key);
        timeInterval(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='NumberChannels';
        Index=strfind(Str,Key);
        numChannels(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        Key='FreqStartInHz';
        Index=strfind(Str,Key);
        FreqStart(i) = sscanf(Str(Index(1) + length(Key):end),'%g',1);
        
        if i > 1
            if ~isequal(numTimeSamples(i),numTimeSamples(1))
                warning('File #%d has a different number of time points than the first file! How are things? Need a hug?\n%s',i,files{i})
            end
            
            if ~isequal(numFrequencies(i),numFrequencies(1))
                warning('File #%d has a different number of frequency points than the first file! Blame it on technology!\n%s',i,files{i})
            end
            
            if ~isequal(startTime(i),startTime(1))
                warning('File #%d has a different start time than the first file! Time''s relative anyways though - amiright?!\n%s',i,files{i})
            end
            
            if ~isequal(timeInterval(i),timeInterval(1))
                warning('File #%d has a different time resolution than the first file! I''m starting to think these data aren''t the data that you think these data are!\n%s',i,files{i})
            end
            
            if ~isequal(FreqInterval(i),FreqInterval(1))
                warning('File #%d has a different frequency resolution than the first file! Science ain''t easy!\n%s',i,files{i})
            end
            
            if ~isequal(numChannels(i),numChannels(1)) %% set extraction based on # of channels
                warning('File #%d has a different number of channels than the first file! BESA hates you AND your timeseries!\n%s',i,files{i})
            end
            
            if ~isequal(numChannels(i),3)
                warning('Ope! File #%d has %.0f channels/orientations! Ya might wanna take a second look at that one! Hint: There should be 3!!!\n%s',i,numChannels(i),files{i})
            end
        end
        
        file_end_freq = end_freq-FreqInterval(i);
        
        
        
        %-----------------------------------------------------------------%
        %                     GET THE DATA                                %
        %-----------------------------------------------------------------%
        j=1;
        
        while ~feof(fid)
            data=fscanf(fid,'%g',numTimeSamples(i));
            
            if ~isempty(data)
                tfc_data(j,:)=data;
                j=j+1;
            end
        end
        [m,~]=size(tfc_data);
        
        %---Close the file---%
        status=fclose(fid);
        
        %-----------------------------------------------------------------%
        %                     AVERAGE TFC DATA                            %
        %-----------------------------------------------------------------%
        
        %---Find the range, in rows, of data to average---%
        freqRowStart = ((start_freq-FreqStart(i))/FreqInterval(i)) + 1;
        freqRowStop = ((file_end_freq-FreqStart(i))/FreqInterval(i)) +1;
        
        %---Average data for average waveforms---%
        for r=1:numChannels(i)
            start_row = numFrequencies(i)*(r-1);
            target_data = tfc_data(start_row+freqRowStart:start_row+freqRowStop,:);
            channel_avg = mean(target_data,1);
            if r == 1
                file_avg = channel_avg;
            else
                file_avg = [file_avg;channel_avg];
            end
        end
        
        
        if ii == 1
            
            all_sub_values(i,1:2,:) = file_avg(2:3,:);
            
        elseif ii == 2
            
            all_sub_values(i,3:4,:) = file_avg(2:3,:);
            
            
        end
        
        
    end
    
end

endTime = startTime(1) + (timeInterval(1)*(size(all_sub_values,3)-1));


time = startTime(1):timeInterval(1):endTime;




answer = questdlg('Would you like to plot the PLV time series for each orientation per subject?', ...
    '', ...
    'Yes','No','Yes');

% Handle response
switch answer
    case 'Yes'
        for i = 1:length(files)
            
            all_ori = squeeze(all_sub_values(i,:,:));
            
            figure(i)
            plot(time,all_ori(1,:))
            hold on
            plot(time,all_ori(2,:))
            hold on
            plot(time,all_ori(3,:))
            hold on
            plot(time,all_ori(4,:))
            
            legend('R11-R21', 'R11-R22','R12-R21', 'R12-R22')
            
            xlabel('Time (ms)')
            ylabel('Amplitude or Power')
            
            ylim([0 Inf]) %This ensures that the y-axis starts at 0
            
        end
    case 'No'
        warning('That is a sad choice...')
end





answer_dom = questdlg('How would you like the dominant orientation to be calculated', ...
    '', ...
    'Average Max Over Time Window','Peak Max Over Time Window','Average Max Over Time Window');


prompt = {'Start Time','End Time'};
title = 'Over what time window should the average or peak be found?';
dims = [1 90];
definput = {'0','0'};
answer_time = inputdlg(prompt,title,dims,definput);
[start_time,end_time] = deal(answer_time{:});
start_time =str2num(start_time);
end_time =str2num(end_time);
end_time = end_time-timeInterval(1);

if start_time >= end_time
    error('enter proper times')
end


%Finding time interval positions

for i = 1:length(time)
    if time(i) == start_time
        initial_time_position = i;
    elseif time(i) == end_time
        final_time_position = i;
    end
end


clear i;
result_table = zeros(n_R11,length(time));
result_table_dominant_or = zeros(n_R11,1);
% Handle response
switch answer_dom
    case 'Average Max Over Time Window'
        for i = 1:length(files)
            for ii = 1:4 %since always 4 options!
                avg_orientation = mean(all_sub_values(i,ii, initial_time_position:final_time_position));
                comparison_vector(ii) = avg_orientation;
                clear avg_orientation;
            end
            if comparison_vector(1) > comparison_vector(2) && comparison_vector(1) > comparison_vector(3) && comparison_vector(1) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,1,:));
                result_table_dominant_or(i,1) = 1;
            elseif comparison_vector(2) > comparison_vector(1) && comparison_vector(2) > comparison_vector(3) && comparison_vector(2) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,2,:));
                result_table_dominant_or(i,1) = 2;
            elseif comparison_vector(3) > comparison_vector(2) && comparison_vector(3) > comparison_vector(1) && comparison_vector(3) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,3,:));
                result_table_dominant_or(i,1) = 3;
            elseif comparison_vector(4) > comparison_vector(2) && comparison_vector(4) > comparison_vector(1) && comparison_vector(4) > comparison_vector(3)
                result_table(i,:) = (all_sub_values(i,4,:));
                result_table_dominant_or(i,1) = 4;
            end
            clear comparison_vector
        end
    case 'Peak Max Over Time Window'
         for i = 1:length(files)
            for ii = 1:4 %since always 4 options!
                max_value_per_orientation = max(all_sub_values(i,ii, initial_time_position:final_time_position));
                comparison_vector(ii) = max_value_per_orientation;
                clear max_value_per_orientation;
            end
            if comparison_vector(1) > comparison_vector(2) && comparison_vector(1) > comparison_vector(3) && comparison_vector(1) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,1,:));
                result_table_dominant_or(i,1) = 1;
            elseif comparison_vector(2) > comparison_vector(1) && comparison_vector(2) > comparison_vector(3) && comparison_vector(2) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,2,:));
                result_table_dominant_or(i,1) = 2;
            elseif comparison_vector(3) > comparison_vector(2) && comparison_vector(3) > comparison_vector(1) && comparison_vector(3) > comparison_vector(4)
                result_table(i,:) = (all_sub_values(i,3,:));
                result_table_dominant_or(i,1) = 3;
            elseif comparison_vector(4) > comparison_vector(2) && comparison_vector(4) > comparison_vector(1) && comparison_vector(4) > comparison_vector(3)
                result_table(i,:) = (all_sub_values(i,4,:));
                result_table_dominant_or(i,1) = 4;
            end
            clear comparison_vector
         end
end


%Finding the dominant orientation based on average or peak

%-----------------------------------------------------------------%
%                     SAVE THE FILE                               %
%-----------------------------------------------------------------%

%---Save the data---%

save_path = uigetdir(path_R11,'Select the directory to the PLV timeseries');
cd(save_path)

%This management of data strcuture sucks... can be highly optimized!
participants = transpose(files);
headers = {'ParID/Dominant Orientation'};
upper_headers = horzcat(cell2table(headers),array2table(time));
output_table = (horzcat(cell2table(participants),array2table(result_table)));

final_table = vertcat(table2cell(upper_headers),table2cell(output_table));
xlswrite('Concatenated_PLV_Dominant_orientations.xlsx',(final_table));

dominant_orientation_table = horzcat(cell2table(participants),array2table(result_table_dominant_or));
xlswrite('Dominant_Orientation_per_File.xlsx',table2cell(dominant_orientation_table));

disp('Done!')