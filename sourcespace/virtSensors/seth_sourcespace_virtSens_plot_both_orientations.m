function seth_sourcespace_virtSens_plot_both_orientations

%PURPOSE:           Plot each orientation for all extracted timeseries.
%                   This allows you to determine if you should vector sum
%                   or pull dominant orientations.
%
%REQUIRED INPUTS:   Fill in the popup window
%		    
%
%NOTES:             
%		   
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, Institute for Human Neuroscience
%VERSION HISTORY:   04/20/2022  v1: First working version of program




%prompt for inputs and convert inputs to number format%
inputs = {'Time Sampling Resolution', 'Epoch Start'};
defaults = {'0', '0'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults, 'on');
[timesample,epochstart] = deal(answer{:});
timesample =str2num(timesample);
epochstart =str2num(epochstart);

%Load files
[files,path,~] = uigetfile('*_Timeseries.tfc','Please Select extracted time series files','Multiselect','on');        %select evt files%

cd(path);


temp = dlmread(files{1},','); %Read one file to calculate time vector


epochend = (epochstart+(timesample*length(temp))-timesample); %Determine ending time


time = linspace(epochstart,epochend,length(temp)); %Create time vector for plotting


for i = 1:length(files)
    
ori_timeseries = dlmread(files{i},',');

ori_timeseries_1 = ori_timeseries(1,:);
ori_timeseries_2 = ori_timeseries(2,:);

figure(i)
plot(time,ori_timeseries_1)
hold on
plot(time,ori_timeseries_2)

legend('Orientation 1', 'Orientation 2')

xlabel('Time (ms)')
ylabel('Amplitude or Power')

ylim([0 Inf]) %This ensures that the y-axis starts at 0

end



end