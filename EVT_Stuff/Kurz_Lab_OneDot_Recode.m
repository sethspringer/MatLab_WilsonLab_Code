function Kurz_Lab_OneDot_Recode

%PURPOSE:           Recode triggers for Kurz lab OneDot PTB data.
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
%31 - 15Hz entrainment onset
%32 - 15Hz entrainment offset
%
%
%RECODED TRIGGERS:
% 101 - 15Hz entrainment onset; TRIAL START ONLY

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
    
    

    
    trial_counter_cond1 = 0;
    trial_counter_cond4 = 0;

    
    %RT array XX
    RT = nan(24,1); %24 possible correct oddball trials
    RT_index = 1;
    
    
    
    
    
    %Loop through the EVT again (stop two from the end, otherwise EVTs ending with fixation codes [12297 or 12298] will have indexing issues)
    for ii = 2:(length(triggers)-2)
        
        %Cond1 - 32Hz
        if (triggers(ii) == (4096+31) || triggers(ii) == (31)) && (triggers(ii-1) == 50 || triggers(ii-1) == (50+4096))
            
            triggers(ii+1) = 101;
            
            trial_counter_cond1 = trial_counter_cond1 + 1;

            %There seems to be an issue where the fixation triggers get dropped rarely. So I need to account for that.
            %The above conditionals should catch most, if there is no
            %fixation, there will be a long gap in time between the last trial ending (32+4096) and the next trial beginning (31)
            
        elseif (triggers(ii) == (4096+31) || triggers(ii) == (31)) && ((time(ii)-time(ii-1)) > 1500000)
            
            triggers(ii+1) = 101;
            
            trial_counter_cond1 = trial_counter_cond1 + 1;
           
            %Cond4 -  Oddball
        elseif (triggers(ii) == (4096+40) || triggers(ii) == (40)) && (triggers(ii-1) == (4096+50) || triggers(ii-1) ==  50 || triggers(ii-1) == (4906 + 32))
            trial_counter_cond4 = trial_counter_cond4+1;
            
            triggers(ii+1) = 102;
            
            %There seems to be an issue where the fixation triggers get dropped rarely. So I need to account for that.
            %The above conditionals should catch most, if there is no
            %fixation, there will be a long gap in time between the last trial ending (32+4096) and the next trial beginning (31)
            
        elseif (triggers(ii) == (4096+40) || triggers(ii) == (40)) && ((time(ii)-time(ii-1)) > 1500000)
            
            triggers(ii+1) = 102;
            
            trial_counter_cond4 = trial_counter_cond4 + 1;
            
 
%             
%              %Cond5 - Response
%         elseif (triggers(ii) == (256) | triggers(ii) == (256+4096)) && (triggers(ii-1) == 4096 | triggers(ii-1) == 256) %left valid cue followed by target%
%             
%             time_since_last_trigger = (time(ii)-time(ii-1))/1000;
%             
%             %Not to correct for some subjects pressing the button during
%             %the entrainment or double pressing, we will do this:
%             
%             
%             if time_since_last_trigger < 100 %there is no way someone can answer in < 100 ms
%                 %Do nothing
%             else
%                 
%                 trial_counter_cond5 = trial_counter_cond5+1;
%                 
%                 RT(RT_index) = time_since_last_trigger;
%                 
%                 RT_index = RT_index + 1;
%             end
%             
%             
        end
        
        
    end
    
    clear ii
%     
%     total_motor = 0;
%     
%      for ii = 2:(length(triggers)-2)
%          if triggers(ii) == 256 | triggers(ii) == (256+4096)
%              total_motor = total_motor+1;
%          end
%      end
%          
     
                       %Writing new EVT file
    evt_info = [time,ones(size(time,1),1),triggers];
    filename = strcat(files{i}(1,1:end-4),'_recoded.evt');
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\n',evt_header);
    fclose(fid);
    dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
         
     
     
     
     
    
    
    
        trial_counter_cond1
        trial_counter_cond4
        

        
      
        
        
        fprintf('\nDone!\n\n')
    
    
end %end of subject for loop

end %end of function
