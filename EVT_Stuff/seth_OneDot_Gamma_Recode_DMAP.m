function seth_OneDot_Gamma_Recode_DMAP

%PURPOSE:           Recode triggers for OneDot Gamma (VisEntrain).
%
%REQUIRED INPUTS:   EVT files for recoding triggers.
%
%		    
%
%NOTES:             Whenever you are writing recodes for these entrainment tasks with a propixx dot, 
%                   there are several trigger sequences that you have to prepare for (example trigger = 20):
%                   1. 20 followed by 20+4096 (the 20+4096 is the propixx dot, lock to this) - MOST COMMON
%                   2. Only 20+4096
%                   3. 20+4096 followed by 20+4096 - UNCOMMON, MOST DANGEROUS; Only lock to the second 20+4096
%		  
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, UNMC
%VERSION HISTORY:   03/01/2022  v2: Fixed fatal flaws with v1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-Recoded Trigger Legend-%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORIGINAL TRIGGERS:
% 31 - 32Hz
% 33 - 40Hz
% 35 - 48Hz
% 37 - Oddball
%
% 50 - Fixation
%
%RECODED TRIGGERS:
% 20 - 32Hz
% 21 - 40Hz
% 22 - 48Hz
% 23 - Oddball
%
% 24 - Fixation

[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%

format long

cd(path);

if ~iscell(files)                                                                        %set number of iterations based on whether there is one or multiple files%
    files = {files};
end
iter = size(files,2);


disp('Processing...');

for i = 1:iter                        %loop throught participants
    %read evt data for file i%
    data = readBESAevt(files{i});
    
    %Separate the triggers and time into vectors that are easier to work on
    triggers = data(:,3);
    time = data(:,1);
    
    %Add not used codes to the data that can be removed later so that the
    %trigger for loop can index past the last position
    triggers(end+1) = 999; 
    time(end+1)     = 999; 

    
    
    
    
    
    trial_counter_cond1 = 0;
    trial_counter_cond2 = 0;
    trial_counter_cond3 = 0;
    trial_counter_cond4 = 0;
    fixation_counter    = 0;
    
    %trial_order = nan(246,2);
    %trial_order_index = 1;
    
    %RT array
    %RT = nan(24,1); %24 possible correct oddball trials
    %RT_index = 1;
    
    
    
    %Loop through the EVT again (stop two from the end, otherwise EVTs ending with fixation codes [12297 or 12298] will have indexing issues)
    for ii = 2:(length(triggers))
        
        %Fixation recode
        if (triggers(ii) == 50+4096 | triggers(ii) == 50) & ((triggers(ii-1) == 50 | triggers(ii-1) == 4096+50 | triggers(ii-1) == 20 | triggers(ii-1) == 21 | triggers(ii-1) == 22 | triggers(ii-1) == 23) & ~(triggers(ii+1) == 4096+50))
            fixation_counter = fixation_counter+1;
            
            triggers(ii) = 24;
            
        %Simple 32Hz recode
        elseif triggers(ii) == 31+4096 & (triggers(ii-1) == 31 | triggers(ii-1) == 4096+31) %Simple triggers
            trial_counter_cond1 = trial_counter_cond1+1;
            
            triggers(ii) = 20;
            
            trial_order(trial_order_index,1) = 1;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1;
                        
            
        %More complicated 32 Hz recode
        elseif triggers(ii) == 4096+31 & (triggers(ii-1) == 24 & ~(triggers(ii+1) == 4096+31)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond1 = trial_counter_cond1+1;
            
            triggers(ii) = 20;
            
            trial_order(trial_order_index,1) = 1;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1;            
            
            
            
        %Simple 40Hz recode
        elseif triggers(ii) == 33+4096 & (triggers(ii-1) == 33 | triggers(ii-1) == 4096+33) %Simple triggers
            trial_counter_cond2 = trial_counter_cond2+1;
            
            
            triggers(ii) = 21;
            
            trial_order(trial_order_index,1) = 2;
            
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1;
                        
            
        %More complicated 40 Hz recode
        elseif triggers(ii) == 4096+33 & (triggers(ii-1) == 24 & ~(triggers(ii+1) == 4096+33)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond2 = trial_counter_cond2+1;
            
            
            triggers(ii) = 21;
            
            trial_order(trial_order_index,1) = 2;
            
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1;       
            
            
            
           
        %Simple 48Hz recode 
        elseif triggers(ii) == 35+4096 & (triggers(ii-1) == 35 | triggers(ii-1) == 4096+35) %Simple triggers
            trial_counter_cond3 = trial_counter_cond3+1;
            
           
            triggers(ii) = 22;
            
            trial_order(trial_order_index,1) = 3;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1; 
                        
            
        %More complicated 48Hz recode
        elseif triggers(ii) == 4096+35 & (triggers(ii-1) == 24 & ~(triggers(ii+1) == 4096+35)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond3 = trial_counter_cond3+1;
            
           
            triggers(ii) = 22;
            
            trial_order(trial_order_index,1) = 3;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            trial_order_index = trial_order_index + 1;                   
            
            
           
            
        %Simple Oddball recode 
        elseif triggers(ii) == 37+4096 & (triggers(ii-1) == 37 | triggers(ii-1) == 4096+37) %Simple triggers
            trial_counter_cond4 = trial_counter_cond4+1;
            
            triggers(ii) = 23;
            
            trial_order(trial_order_index,1) = 4;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            
            trial_order_index = trial_order_index + 1;    
                        
            
        %More complicated Oddball recode
        elseif triggers(ii) == 4096+37 & ((triggers(ii-1) == 24 | triggers(ii-1) == 256) & ~(triggers(ii+1) == 4096+37)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond4 = trial_counter_cond4+1;
            
            triggers(ii) = 23;
            
            trial_order(trial_order_index,1) = 4;
            trial_order(trial_order_index,2) = time(ii)/1000000;
            
            
            trial_order_index = trial_order_index + 1;    
            
            
            
            
        elseif triggers(ii) == 256 & triggers(ii-1) == 24
            
            RT(RT_index) = (time(ii)-time(ii-1))/1000;
            
            RT_index = RT_index + 1;
            

            
        end %End of IFs/ELSEIFs
        
        
    end %End of trigger for loop
    
    
    %Delete the 999s that were added for padding
    time(triggers==999) = [];
    triggers(triggers==999) = [];
    
    
    
        trial_counter_cond1
        trial_counter_cond2
        trial_counter_cond3
        trial_counter_cond4
        
        fixation_counter
        
        trial_order
        RT
        
        ave_RT = mean(RT)
        
                  %Writing new EVT file
    evt_info = [time,ones(size(time,1),1),triggers];
    filename = strcat(files{i}(1,1:end-4),'_recoded.evt');
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\n',evt_header);
    fclose(fid);
    dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
        
        
        fprintf('\nDone with %s!\n\n',files{i}(1,1:3))
    
    
end
