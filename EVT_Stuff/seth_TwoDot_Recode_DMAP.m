function seth_TwoDot_Recode_DMAP

%PURPOSE:           Recode triggers for DMAP TwoDot.
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
%AUTHOR:            Seth D. Springer, DICoN Lab, Boys Town National Research Hospital
%VERSION HISTORY:   03/01/2022  v2: Fixed fatal flaws with v1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-Recoded Trigger Legend-%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORIGINAL TRIGGERS:
% 21 - Attend Left
% 23 - Left Oddball
% 25 - Attend Right
% 27 - Right Oddball
%
% 50 - Fixation
%
%RECODED TRIGGERS:
% 31 - Attend Left
% 32 - Left Oddball
% 33 - Attend Right
% 34 - Right Oddball
%
% 24 - Fixation


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%



cd(path);

if ~iscell(files)                                                                        %set number of iterations based on whether there is one or multiple files%
    files = {files};
end
iter = size(files,2);


output_table = [];
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
    trial_counter_cond5 = 0;
    
    fixation_counter = 0;
    
    
    %RT array
    RT = nan(12,1); %12 possible correct oddball trials
    RT_index = 1;
    
    
    
    %Loop through the EVT again (stop two from the end, otherwise EVTs ending with fixation codes [12297 or 12298] will have indexing issues)
    for ii = 2:(length(triggers))
        
        
        %Fixation recode
        if (triggers(ii) == 50+4096 | triggers(ii) == 50) & ((triggers(ii-1) == 50 | triggers(ii-1) == 4096+50 | triggers(ii-1) == 31 | triggers(ii-1) == 32 | triggers(ii-1) == 33 | triggers(ii-1) == 34) & ~(triggers(ii+1) == 4096+50))
            fixation_counter = fixation_counter+1;
            
            triggers(ii) = 24;
            
            
        %Simple Attend Left recode
        elseif triggers(ii) == 21+4096 & (triggers(ii-1) == 21 | triggers(ii-1) == 4096+21) %Simple triggers
            trial_counter_cond1 = trial_counter_cond1+1;
            
            triggers(ii) = 31;
            
                        
            
        %More complicated Attend Left recode
        elseif triggers(ii) == 4096+21 & (triggers(ii-1) == 24 & ~(triggers(ii+1) == 4096+21)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond1 = trial_counter_cond1+1;
            
            triggers(ii) = 31;
                   
            
            
        %Simple Left Oddball recode
        elseif triggers(ii) == 23+4096 & (triggers(ii-1) == 23 | triggers(ii-1) == 4096+23) %Simple triggers
            trial_counter_cond2 = trial_counter_cond2+1;
            
            
            triggers(ii) = 32;
            
                        
            
        %More complicated Left Oddball recode
        elseif triggers(ii) == 4096+23 & ((triggers(ii-1) == 24 | triggers(ii-1) == 256  | triggers(ii-1) == 512) & ~(triggers(ii+1) == 4096+23)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond2 = trial_counter_cond2+1;
            
            
            triggers(ii) = 32;
            
           
            
        %Simple Attend Right recode 
        elseif triggers(ii) == 25+4096 & (triggers(ii-1) == 25 | triggers(ii-1) == 4096+25) %Simple triggers
            trial_counter_cond3 = trial_counter_cond3+1;
            
           
            triggers(ii) = 33;
            
                        
            
        %More complicated Attend Right recode
        elseif triggers(ii) == 4096+25 & (triggers(ii-1) == 24 & ~(triggers(ii+1) == 4096+25)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond3 = trial_counter_cond3+1;
            
           
            triggers(ii) = 33;
            
            
        %Simple Right Oddball recode 
        elseif triggers(ii) == 27+4096 & (triggers(ii-1) == 27 | triggers(ii-1) == 4096+27) %Simple triggers
            trial_counter_cond4 = trial_counter_cond4+1;
            
            triggers(ii) = 34;
            
                        
            
        %More complicated Right Oddball recode
        elseif triggers(ii) == 4096+27 & ((triggers(ii-1) == 24 | triggers(ii-1) == 256 | triggers(ii-1) == 512) & ~(triggers(ii+1) == 4096+27)) %Ensures that the transformed triggers (x+4096) only count the second one, but also are counted if there is only one trigger and not the regular pair
            trial_counter_cond4 = trial_counter_cond4+1;
            
            triggers(ii) = 34;

            
        elseif (triggers(ii) == 256 | triggers(ii) == 512) & (triggers(ii-1) == 34 | triggers(ii-1) == 32)
            
            RT(RT_index) = (time(ii)-time(ii-1))/1000;
            
            RT_index = RT_index + 1;        
            
        end
        
        
    end
    

    %Delete the 999s that were added for padding
    time(triggers==999) = [];
    triggers(triggers==999) = [];
    
    
    
         
     
                       %Writing new EVT file
    evt_info = [time,ones(size(time,1),1),triggers];
    filename = strcat(files{i}(1,1:end-4),'_recoded.evt');
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\n',evt_header);
    fclose(fid);
    dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
         
     
     
     
     
    
    
    
        trial_counter_cond1
        trial_counter_cond2
        trial_counter_cond3
        trial_counter_cond4
        
        
        fixation_counter
        
        RT

        
      
        
        
        fprintf('\nDone!\n\n')
    
    
end %end of subject for loop

end %end of function
