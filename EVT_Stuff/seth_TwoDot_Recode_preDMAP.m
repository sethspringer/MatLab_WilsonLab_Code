function seth_TwoDot_Recode_preDMAP

%PURPOSE:           Recode triggers for Pre-DMAP TwoDot (PESCAH and SCAN Pilot).
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
%
%RECODED TRIGGERS:
% 31 - Attend Left
% 32 - Left Oddball
% 33 - Attend Right
% 34 - Right Oddball
%


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
    
    
    
    %Delete nonsense triggers 12303 and 12302 from data%
    %time(triggers==4096) = [];
    %triggers(triggers==4096) = [];
    
    
    
    
    trial_counter_cond1 = 0;
    trial_counter_cond2 = 0;
    trial_counter_cond3 = 0;
    trial_counter_cond4 = 0;
    trial_counter_cond5 = 0;
    
    %RT array
    RT = nan(24,1); %24 possible correct oddball trials
    RT_index = 1;
    
    
    
    
    
    %Loop through the EVT again (stop two from the end, otherwise EVTs ending with fixation codes [12297 or 12298] will have indexing issues)
    for ii = 2:(length(triggers)-2)
        
        %Cond1 - Left Attend
        if triggers(ii) == (4096+21) && (triggers(ii-1) == 4096 | triggers(ii-1) ==  21 | triggers(ii-1) == 256 | triggers(ii-1) == (4096+21)) %left valid cue followed by target%
            trial_counter_cond1 = trial_counter_cond1+1;
            
            %Now we have to check to make sure that the following trigger isn't the propixx trigger, otherwise, we want to lock to that...
            if triggers(ii+1) == (4096+21)
                triggers(ii+1) = 31;
            else
                triggers(ii) = 31;
            end
            

            %Cond2 -  Left Oddball
        elseif triggers(ii) == (4096+23) && (triggers(ii-1) == 4096 | triggers(ii-1) ==  23 | triggers(ii-1) == 256 | triggers(ii-1) == (4096+23)) %left valid cue followed by target%
            trial_counter_cond2 = trial_counter_cond2+1;
                        
            %Now we have to check to make sure that the following trigger isn't the propixx trigger, otherwise, we want to lock to that...
            if triggers(ii+1) == (4096+23)
                triggers(ii+1) = 32;
            else
                triggers(ii) = 32;
            end            
            
            %Cond3 - Right Attend
        elseif triggers(ii) == (4096+25) && (triggers(ii-1) == 4096 | triggers(ii-1) ==  25 | triggers(ii-1) == 256 | triggers(ii-1) == (4096+25)) %left valid cue followed by target%
            trial_counter_cond3 = trial_counter_cond3+1;
            
            %Now we have to check to make sure that the following trigger isn't the propixx trigger, otherwise, we want to lock to that...
            if triggers(ii+1) == (4096+25)
                triggers(ii+1) = 33;
            else
                triggers(ii) = 33;
            end            
            
            
           
            %Cond4 -  Right Oddball
        elseif triggers(ii) == (4096+27) && (triggers(ii-1) == 4096 | triggers(ii-1) ==  27 | triggers(ii-1) == 256 | triggers(ii-1) == (4096+27)) %left valid cue followed by target%
            trial_counter_cond4 = trial_counter_cond4+1;
            
            %Now we have to check to make sure that the following trigger isn't the propixx trigger, otherwise, we want to lock to that...
            if triggers(ii+1) == (4096+27)
                triggers(ii+1) = 34;
            else
                triggers(ii) = 34;
            end            
            
             %Cond5 - Response
        elseif (triggers(ii) == (256) | triggers(ii) == (256+4096)) && (triggers(ii-1) == 4096 | triggers(ii-1) == 256) %left valid cue followed by target%
            trial_counter_cond5 = trial_counter_cond5+1;
            
            
            %Need to find the last oddball flicker start
            all_oddball_locations = find(triggers == 32 | triggers == 34);
            
            most_recent_oddball_start = all_oddball_locations(end);
            
            
            RT(RT_index) = (time(ii)-time(most_recent_oddball_start))/1000;
            
            RT_index = RT_index + 1;    
            
            
        end
        
        
    end
    
    clear ii
    
    total_motor = 0;
    
     for ii = 2:(length(triggers)-2)
         if triggers(ii) == 256 | triggers(ii) == (256+4096)
             total_motor = total_motor+1;
         end
     end
         
     
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
        trial_counter_cond5
        total_motor
        
        RT

        
      
        
        
        fprintf('\nDone!\n\n')
    
    
end %end of subject for loop

end %end of function
