function sensorspace_plot_one_channel_plot_all(start_time, end_time, start_freq, end_freq, low_limit_scale, high_limit_scale) 

%PURPOSE:           Plot multiple TFCs (use after sensorspace_plot_one_channel)
%
%REQUIRED INPUTS:   Single sensor TFC data.
%
%
%
%NOTES:             YOU MUST have saved ColorScheme information using the
%                   information on the wiki.
%                   For an example of this look at D:\MIND_VWM\derivative\SensorSpace\1Hz50ms__4-50hz\RELpwr\GenMTG
%                   on my computer (BT137868)
%
%
%AUTHOR:            Seth D. Springer, IHN, University of Nebraska Medical Center
%VERSION HISTORY:   09/24/2021  v1: First working version of program
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[files,path,~] = uigetfile('*.tfc','Please Select TFCs','Multiselect','on');        %select evt files%

cd(path)

load('ColorScheme') %Load in the ColorScheme.mat that you made and want to apply to all figures


if ~iscell(files)             
    files = {files};
end


for i = 1:length(files)
    
    %Grab naming stuff to edit the figure name and titles
    filename    = files{i};
    sensor_name = filename(1:7);
    parID       = filename(9:11);
    
    figure; %Open new figure for each subject
    
    besa_tftplot(filename,sensor_name,[start_time end_time], [start_freq end_freq]);
    
    set(gcf, 'Name',sensor_name);
    
    %Automatically change axis and color scheme based on a color scheme that you created before
    set(gca, 'Colormap', variable_name, 'CLim', [low_limit_scale high_limit_scale]); %'variable_name' was created when you made the ColorScheme
    
    title(parID)
    
end

end
