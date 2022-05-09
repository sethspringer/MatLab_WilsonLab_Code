function Seth_Extract_TFC_TimeseriesV3(startFreq,endFreq)

%Normal Extract_TFC... function but output files are moved into a new
%"Timeseries" folder at the end.

%---Load All Files---%
files = dir('*.tfc');

for i=1:length(files)
    [~,name,ext] = fileparts(files(i).name);
    filename = strcat(name,ext);
    newname = strcat(name,'_TFC_Timeseries',ext);

%---Create Average Waveform For Each File---%
    Str=fileread(filename);
    fid=fopen(filename,'rt');

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
                warning('File #%d has a different number of time points than the first file! How are things? Need a hug?\n%s',i,files(i).name)
            end
            
            if ~isequal(numFrequencies(i),numFrequencies(1))
                warning('File #%d has a different number of frequency points than the first file! Blame it on technology!\n%s',i,files(i).name)
            end
            
            if ~isequal(startTime(i),startTime(1))
                warning('File #%d has a different start time than the first file! Time''s relative anyways though - amiright?!\n%s',i,files(i).name)
            end
            
            if ~isequal(timeInterval(i),timeInterval(1))
                warning('File #%d has a different time resolution than the first file! I''m starting to think these data aren''t the data that you think these data are!\n%s',i,files(i).name)
            end
            
            if ~isequal(FreqInterval(i),FreqInterval(1))
                warning('File #%d has a different frequency resolution than the first file! Science ain''t easy!\n%s',i,files(i).name)
            end
            
            if ~isequal(numChannels(i),numChannels(1)) %% set extraction based on # of channels
                warning('File #%d has a different number of channels than the first file! BESA hates you AND your timeseries!\n%s',i,files(i).name)
            end
            
            if ~isequal(numChannels(i),2)
                warning('Ope! File #%d has %.0f channels/orientations! Ya might wanna take a second look at that one!\n%s',i,numChannels(i),files(i).name)
            end
        end
        
        file_endFreq = endFreq-FreqInterval(i);


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
    [m,n]=size(tfc_data);
    
    %---Close the file---%
    status=fclose(fid);

        %-----------------------------------------------------------------%
        %                     AVERAGE TFC DATA                            %
        %-----------------------------------------------------------------%

        %---Find the range, in rows, of data to average---%
        freqRowStart = ((startFreq-FreqStart(i))/FreqInterval(i)) + 1;
        freqRowStop = ((file_endFreq-FreqStart(i))/FreqInterval(i)) +1;

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
    %     figure(1)
    %     plot(avgWave(1,:))
    %     figure(2)
    %     plot(tfc_data(1,:))
        %-----------------------------------------------------------------%
        %                     SAVE THE FILE                               %
        %-----------------------------------------------------------------%

        %---Save the data---%
        disp('Saving data...')

        csvwrite(newname,file_avg);
        
        clearvars -except files startFreq endFreq numTimeSamples FreqInterval numFrequencies startTime timeInterval numChannels FreqStart
    
end

movefile('*Timeseries*','Timeseries')
movefile('

disp('Done!')



