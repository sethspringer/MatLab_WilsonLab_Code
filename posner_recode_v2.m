function [behavior] = posner_recode_v2(varargin)

%PURPOSE:           Recode triggers for posner.
%
%REQUIRED INPUTS:   EVT files for recoding MAMMOA posner triggers - updated.
%
%		    
%
%NOTES:            Removes ALL 12303/12302 triggers
%                  Removes any 12299 triggers that are not preceded by a
%                  response (128** or 125**) or by 12289-12296 (all targets)
%                  Removes any 12295 triggers that are not preceded by
%                  12297 or 12298
%		   
%
%
%                  
%AUTHOR:            Alex I. Wiesman, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   10/30/2018  v1: First working version of program
%                   04/03/2020  v2: Now removes all misc triggers as stated
%                   above, but also recodes target triggers (see legend
%                   below) - also now outputs behavior
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-Recoded Trigger Legend-%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15 = Correct Valid Cue
% 16 = Correct Invalid Cue
% 25 = Correct Valid Target
% 26 = Correct Invalid Target
% 35 = Correct Valid Response
% 36 = Correct Invalid Response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 12297 or 12298 = cue
% 12289 - 12296 are targets (89, 90, 93, 94 -> 128*; 91, 92, 95, 96 -> 125*)
% 128* and 125* are then responses
% valid targets are 12289, 94, 91, 96
% invalid targets are 12290, 93, 92, 95

if isempty(varargin)
    sd_cutoff = 1000;
else
    sd_cutoff = varargin{1};
end

[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%
    
cd(path);  





if isa(files,'char')                                                                        %set number of iterations based on whether there is one or multiple files%
    iter = 1;
else
    iter = length(files);
end



%Initialization
behavior.group_ACC = NaN(1,iter);
behavior.group_valid_ACC = NaN(1,iter);
behavior.group_invalid_ACC = NaN(1,iter);

behavior.group_RT = NaN(1,iter);
behavior.group_valid_RT = NaN(1,iter);
behavior.group_invalid_RT = NaN(1,iter);

behavior.group_RT_OE = NaN(1,iter);
behavior.group_valid_RT_OE = NaN(1,iter);
behavior.group_invalid_RT_OE = NaN(1,iter);



for i = 1:iter                        %loop throught participants                  
    if isa(files,'char')                                                                    %read evt data for file i%
        data = readBESAevt(files);
    elseif length(files) > 1;
        data = readBESAevt(files{i});
    else
        error('No compatible files selected! Boom roasted!')
    end
    
    
    
    
    %Separate the triggers and time into vectors that are easier to work on
    triggers = data(:,3);
    time = data(:,1);
    
    
    
    %Delete nonsense triggers 12303 and 12302 from data%
    time(triggers==12303) = [];                     
    triggers(triggers==12303) = [];
    time(triggers==12302) = [];                     
    triggers(triggers==12302) = [];
    time(triggers==0) = [];                     
    triggers(triggers==0) = [];
    
    
  
    %Loop through the EVT
    for ii = 3:length(triggers)
        if triggers(ii) == 12299
            prev_trig = num2str(triggers(ii-1));
            
            
            if strcmp(prev_trig(1:3),'125') | strcmp(prev_trig(1:3),'128') | strcmp(prev_trig,'12289') | strcmp(prev_trig,'12290') | strcmp(prev_trig,'12291') | strcmp(prev_trig,'12292') | strcmp(prev_trig,'12293') | strcmp(prev_trig,'12294') | strcmp(prev_trig,'12295') | strcmp(prev_trig,'12296')
                                                    %keep these%
            else
                triggers(ii) = NaN;
                time(ii) = NaN;
            end
        elseif triggers(ii) == 12295
            prev_trig = num2str(triggers(ii-1));
            if strcmp(prev_trig,'12297') | strcmp(prev_trig,'12298')
                %keep these%
            else
                triggers(ii) = NaN;
                time(ii) = NaN;
            end
        end
    end
    triggers(isnan(triggers)) = [];
    time(isnan(time)) = [];
    
    
    
    
    
    
    
    valid_trial_counter = 0;
    invalid_trial_counter = 0;
    RT = NaN(1,length(triggers));
    validity_index = NaN(1,length(triggers));
    for ii = 3:length(triggers)
        if (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12289 | triggers(ii+1) == 12294 | triggers(ii+1) ==  12291 | triggers(ii+1) == 12296) %valid cue followed by target%
            valid_trial_counter = valid_trial_counter+1;
            if ((triggers(ii+1) == 12289 | triggers(ii+1) == 12294) && (triggers(ii+2) > 12799 && triggers(ii+2) < 12900)) | ((triggers(ii+1) == 12291 | triggers(ii+1) == 12296) && (triggers(ii+2) > 12499 && triggers(ii+2) < 12600)) %correct valid target%
                triggers(ii) = 15;
                triggers(ii+1) = 25;
                triggers(ii+2) = 35;
                validity_index(1,ii) = 1;
                RT(1,ii) = ((time(ii+2)-time(ii+1))-48000)/1000;
            end
        elseif (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12290 | triggers(ii+1) == 12293 | triggers(ii+1) ==  12292 | triggers(ii+1) == 12295) %invalid cue followed by target%
            invalid_trial_counter = invalid_trial_counter+1;
            if ((triggers(ii+1) == 12290 | triggers(ii+1) == 12293) && (triggers(ii+2) > 12799 && triggers(ii+2) < 12900)) | ((triggers(ii+1) == 12292 | triggers(ii+1) == 12295) && (triggers(ii+2) > 12499 && triggers(ii+2) < 12600)) %correct invalid target%                if triggers(ii+1) == 12291 | triggers(ii+1) == 12296 %correct valid trial%
                triggers(ii) = 16;
                triggers(ii+1) = 26;
                triggers(ii+2) = 36;
                validity_index(1,ii) = 0;
                RT(1,ii) = ((time(ii+2)-time(ii+1))-48000)/1000;
            end
        end
    end
    
    RT(isnan(RT)) = [];
    validity_index(isnan(validity_index)) = []; 
    
    if isa(files,'char')
        fprintf('Recoded %d trials for file %s\n',size(RT,2),files);
    else
        fprintf('Recoded %d trials for file %s\n',size(RT,2),files{i});
    end
    
    if size(RT,2) > 1 && size(validity_index,2) > 1
        behavior.group_ACC(1,i) = size(validity_index,2)/(valid_trial_counter+invalid_trial_counter);
        behavior.group_valid_ACC(1,i) = sum(validity_index == 1)/valid_trial_counter;
        behavior.group_invalid_ACC(1,i) = sum(validity_index == 0)/invalid_trial_counter;
        
        validity_index = logical(validity_index);
        behavior.group_RT(1,i) = mean(RT,2);
        behavior.group_valid_RT(1,i) = mean(RT(validity_index),2);
        behavior.group_invalid_RT(1,i) = mean(RT(~validity_index),2);
        
        RT_OE = RT(RT<(mean(RT)+(sd_cutoff*std(RT))) & RT>(mean(RT)-(sd_cutoff*std(RT))));
        validity_index_OE = validity_index(RT<(mean(RT)+(sd_cutoff*std(RT))) & RT>(mean(RT)-(sd_cutoff*std(RT))));
        behavior.group_RT_OE(1,i) = mean(RT_OE,2);
        behavior.group_valid_RT_OE(1,i) = mean(RT_OE(validity_index_OE),2);
        behavior.group_invalid_RT_OE(1,i) = mean(RT_OE(~validity_index_OE),2);
    end
    

        evt_info = [time,data(1:size(time,1),2),triggers];
        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_recoded.evt');
        else
            filename = strcat(files{i}(1,1:end-4),'_recoded.evt');                              %save evt file%
        end
        fid = fopen(filename,'wt');
        fprintf(fid,'%s\n',evt_header);
        fclose(fid);
        dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
end