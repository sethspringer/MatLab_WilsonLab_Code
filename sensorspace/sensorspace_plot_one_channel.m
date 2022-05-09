function sensorspace_plot_one_channel(channel_to_plot)

%PURPOSE:           Grab and save data for one sensor of a TFC for each person in a file list.
%
%REQUIRED INPUTS:   GenMTG REL TFC (can do ABS but wouldn't be useful usually)
%
%                   You can then use sensorspace_plot_one_channel_plot_all
%                   to see the plots for the sensor or interest per person
%
%NOTES:             The channel that you enter needs the "MEG" prefix
%
%
%AUTHOR:            Seth D. Springer, IHN, University of Nebraska Medical Center
%VERSION HISTORY:   09/24/2021  v1: First working version of program
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[files,path,~] = uigetfile('*.tfc','Please Select TFCs','Multiselect','on');        %select evt files%

cd(path)

% defaults
GA.ConditionName = '';
GA.DataType = 'ERDERS_AMP';
GA.NumberOfTrials = 1;
GA.StatisticsCorrection='off';
GA.EvokedSignalSubtraction='off';


%# iterations
n     = length(files);


% read first file to get info on channels, time and frequency bins
hdr   = readBESAtfc(files{1});
nchan = size(hdr.Data,1);
nfreq = length(hdr.Frequency);
ntime = length(hdr.Time);


% read data into array
dataArr = zeros(n,nchan,ntime,nfreq);
for ii=1:n
    tmp = readBESAtfc(files{ii});
    try
        dataArr(ii,:,:,:)   = tmp.Data;
        fprintf('Grabbing data from file # %.0f of %.0f\n',ii,length(files));
    catch errorcode
        fprintf('Data mismatch error on file # %.0f \n',ii);
    end
end

%In 'dataArr':
%position 1 = # of files             (80)
%position 2 = # of sensors           (173)
%position 3 = # of time measurements (289)
%position 4 = # of freq measurements (24)


%Initalizing counter to make sure only 1 channel is found
count = 0;

%Loop throught to find the channel position and get rid of other channels
for iii=1:nchan
    
    %Grabbing channel names for the search
    GA.ChannelLabels = hdr.ChannelLabels(iii,:);
    GA.ConditionName = strtrim(GA.ChannelLabels);
    
    if strcmp(GA.ConditionName,channel_to_plot)   %Look for the channel of interest only
        
        dataArr_oneSensor = dataArr(:,iii,:,:);
        
        channel_index = iii;
        
        count=count+1;
        
    end
end

%Check if more than one channel was found, otherwise, move on.
if count > 1
    error('The function found too many channels, you show look into this...')
end


dataArr_oneSensor = squeeze(dataArr_oneSensor); %Remove the "sensor" dimension bc there is only one now

%In 'dataArr_oneSensor':
%position 1 = # of files             (80)
%position 2 = # of time measurements (289)
%position 3 = # of freq measurements (24)


%Initialize
GA.Data = zeros(1,ntime,nfreq);

%Setting the sensor name, one last time
GA.ChannelLabels = hdr.ChannelLabels(channel_index,:);
GA.ConditionName = strtrim(GA.ChannelLabels);



% Loop through each person to save out their channel of interest
for iiii=1:n
    
    subjectID = (files{iiii}(1:3));
    GA.ConditionName = strtrim(GA.ChannelLabels); %Reseting the name so that it is not elongated
    GA.ConditionName = strcat(GA.ConditionName,'_',subjectID);
    
    
    GA.Data(1,:,:) = dataArr_oneSensor(iiii,:,:);  %Assign the iiiith subject's data
    GA.Time = hdr.Time;
    GA.Frequency = hdr.Frequency;
    besa_writetfc(GA);
    
end


%Create average for comparison
GA.ConditionName = strtrim(GA.ChannelLabels); %Reseting the name
GA.ConditionName = strcat(GA.ConditionName,'_average');

dataArr_oneSensor_average = squeeze(mean(dataArr_oneSensor,1));

GA.Data(1,:,:) = dataArr_oneSensor_average(:,:);  %Assign the iiiith subject's data
GA.Time = hdr.Time;
GA.Frequency = hdr.Frequency;
besa_writetfc(GA);


end
