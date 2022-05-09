function extractPeakFreq(FrequencyRange,FrequencyResolution,LowestFreq)

%PURPOSE:           Extract the peak freq on the txt file that you got from
%                   running the "sourcespace_virtSens_TFA_extract_frequencySeries"function
%
%REQUIRED INPUTS:   Txt file made my the above mentioned function.
%                   FrequencyRange variable in the form of [X Y], in which
%                   X is the lowest frequency of interest and Y is the
%                   highest frequency of interest.
%                   EXAMPLE: Want to exact peak alpha (7-12 Hz); [7 12]
%                   
%                   Also need to input what is the frequency resolution and
%                   lowest frequency of the freqSeries that is input.
%		    
%
%NOTES:             HALFWAY THROUGH WRITING THIS I REALIZED THAT THIS
%                   ALREADY EXISTS: sourcespace_virtSens_TFC_find_peakFreq
%		  
%                  
%AUTHOR:            Seth D Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   11/9/2021  v1: First working version of program


%%%%   Testing    %%%%

FrequencyRange = [4 14];

FrequencyResolution = .5;

LowestFreq = 4;

%%%%%%%%%%%%%%%%%%%%%%


[file,path,~] = uigetfile('*.txt','Please Select The FreqSeries Text File');   

filepath = fullfile(path,file);
opt = detectImportOptions(filepath);

data = readtable(filepath,opt);
parID = table2cell(data(:,1));
values = table2array(data(:,2:end));

num_of_values   = size(values,2);
num_of_subjects = size(values,1);

%Need to creat a vector to match input values with their frequencies
frequencies = [LowestFreq:FrequencyResolution:200];  %make it too long
actual_freqs = frequencies(1:num_of_values); %Then trim it up

%Index the positions where the start and end frequencies are
freq_range_of_interest(1) = find(actual_freqs==FrequencyRange(1));
freq_range_of_interest(2) = find(actual_freqs==FrequencyRange(2));


for i = 1:num_of_subjects
    
    loop_values = values(i,:);
    
    loop_values_for_freq_of_interest = loop_values(freq_range_of_interest(1):freq_range_of_interest(2));
    
    
    
    
end













































end