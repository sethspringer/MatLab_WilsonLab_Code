function posner_recode_v3

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
%AUTHOR:            Alex I. Wiesman, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   10/30/2018  v1: First working version of program
%                   04/03/2020  v2: Now removes all misc triggers as stated
%                   above, but also recodes target triggers (see legend
%                   below) - also now outputs behavior
%                   
%                   Seth D. Springer, DICoN Lab, UNMC
%                   06/29/2021  v3: Updated with more comments and added
%                   separate triggers for left vs. right targets for
%                   saccade analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-Recoded Trigger Legend-%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15 = Correct Left Valid Cue
% 16 = Correct Left Invalid Cue
% 25 = Correct Left Valid Target
% 26 = Correct Left Invalid Target
% 35 = Correct Left Valid Response
% 36 = Correct Left Invalid Response

% 18 = Correct Right Valid Cue
% 19 = Correct Right Invalid Cue
% 28 = Correct Right Valid Target
% 29 = Correct Right Invalid Target
% 38 = Correct Right Valid Response
% 39 = Correct Right Invalid Response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 12297 or 12298 = cue
% 12289 - 12296 are targets (89, 90, 93, 94 -> 128*; 91, 92, 95, 96 -> 125*)
% 128* and 125* are then responses
% valid left  targets are 12289 and 91
% valid right targets are 12294 and 96
%
% invalid left  targets are 12293 and 95
% invalid right targets are 12290 and 92


% if isempty(varargin)
%     sd_cutoff = 1000;
% else
%     sd_cutoff = varargin{1};
% end

inputs = {'SD Cutoff'};
defaults = {'2.5'};
answer = inputdlg(inputs, 'What SD Cutoff should be used for RT outlier exclusions', 2, defaults, 'on');
sd_cutoff = str2double(answer);


if isnan(sd_cutoff)
    error("You suck, enter a number...")
end


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%

if isempty(files)
    error('No compatible files selected! Boom roasted!')
end


cd(path);  

if ~iscell(files)                                                                        %set number of iterations based on whether there is one or multiple files%
    files = {files};
end
iter = size(files,2);


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



    %%%%%%%%%%%%%%%%%%% Creating .xls with behavioral data %%%%%%%%%%%%%%%%%%%%%
DefaultName = 'Posner_Recode_Behavior_Data';
[FileName,PathName,~] = uiputfile('*.csv','Please select path for output behavioral file',DefaultName);
cd(PathName);
%outputfile = fopen(FileName,'w+'); 						%open a file for writing participant data and write headers%
%fprintf(outputfile, 'ParID\t total_correct_trial_number\t group_ACC\t group_valid_ACC\t group_invalid_ACC\t group_RT\t group_valid_RT\t group_invalid_RT\t group_RT_OE\t group_valid_RT_OE\t group_invalid_RT_OE\n');        %Headers applied to the output file


output_table = [];
disp('Processing...');
for i = 1:iter                        %loop throught participants                  
                                                                  %read evt data for file i%
    data = readBESAevt(files{i});
    
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
            
             %keep 12299s that are preceeded by these codes%
            if strcmp(prev_trig(1:3),'125') | strcmp(prev_trig(1:3),'128') | strcmp(prev_trig,'12289') | strcmp(prev_trig,'12290') | strcmp(prev_trig,'12291') | strcmp(prev_trig,'12292') | strcmp(prev_trig,'12293') | strcmp(prev_trig,'12294') | strcmp(prev_trig,'12295') | strcmp(prev_trig,'12296')
                      
            else
                %delete 12299s that dont follow the above codes%
                triggers(ii) = NaN;
                time(ii) = NaN;
            end
            
            %Remove 12295s that arent preceeded by cues (12297/12298)
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
    
    %Remove marked positions from the For loop above%
    triggers(isnan(triggers)) = [];
    time(isnan(time)) = [];
    
    
    
    
    
    
    %Initialize the trial counters%
    valid_trial_counter = 0;
    invalid_trial_counter = 0;
    
    %Initialize the RT and validity vectors%
    RT = NaN(1,length(triggers));
    validity_index = NaN(1,length(triggers));
    
    
    %Loop through the EVT again (stop two from the end, otherwise EVTs ending with fixation codes [12297 or 12298] will have indexing issues)
    for ii = 3:(length(triggers)-2)
        
        %Left Valid
        if (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12289 | triggers(ii+1) ==  12291) %left valid cue followed by target%
            valid_trial_counter = valid_trial_counter+1;
            
            if (triggers(ii+1) == 12289 & triggers(ii+2) == 12801) | (triggers(ii+1) == 12291 & triggers(ii+2) == 12547)
            
            triggers(ii)   = 15;
            triggers(ii+1) = 25;
            triggers(ii+2) = 35;
            
            validity_index(1,ii) = 1;
            %4800ms is the screen delay
            RT(1,ii) = ((time(ii+2)-time(ii))-48000)/1000;
            end
            
        %Right Valid    
        elseif (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12294 | triggers(ii+1) ==  12296) %right valid cue followed by target%
            valid_trial_counter = valid_trial_counter+1;
            
            if (triggers(ii+1) == 12294 & triggers(ii+2) == 12806) | (triggers(ii+1) == 12296 & triggers(ii+2) == 12552)
            
            triggers(ii)   = 18;
            triggers(ii+1) = 28;
            triggers(ii+2) = 38;
            
            validity_index(1,ii) = 1;
            RT(1,ii) = ((time(ii+2)-time(ii))-48000)/1000;
            end
           
        %Left Invalid
        elseif (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12293 | triggers(ii+1) ==  12295) %right valid cue followed by target%
            invalid_trial_counter = invalid_trial_counter+1;
            if (triggers(ii+1) == 12293 & triggers(ii+2) == 12805) | (triggers(ii+1) == 12295 & triggers(ii+2) == 12551)
                
                triggers(ii)   = 16;
                triggers(ii+1) = 26;
                triggers(ii+2) = 36;
                
                validity_index(1,ii) = 0;
                RT(1,ii) = ((time(ii+2)-time(ii))-48000)/1000;
            end
            
            
            
            
        %Right Invalid
        elseif (triggers(ii) == 12297 | triggers(ii) == 12298) && (triggers(ii+1) == 12290 | triggers(ii+1) ==  12292) %right valid cue followed by target%
            invalid_trial_counter = invalid_trial_counter+1;
            if (triggers(ii+1) == 12290 & triggers(ii+2) == 12802) | (triggers(ii+1) == 12292 & triggers(ii+2) == 12548)
                
                triggers(ii)   = 19;
                triggers(ii+1) = 29;
                triggers(ii+2) = 39;
                
                validity_index(1,ii) = 0;
                RT(1,ii) = ((time(ii+2)-time(ii))-48000)/1000;
            end
            
            
        end
    end
    
    
    RT(isnan(RT)) = [];
    validity_index(isnan(validity_index)) = [];
    
    
    %Behavioral calculations
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
    

    
   %Print the behavioral data to an excel sheet

   headers = {'ParID' 'total_correct_trial_number' 'group_ACC' 'group_valid_ACC' 'group_invalid_ACC' 'group_RT'...
              'group_valid_RT' 'group_invalid_RT' 'group_RT_OE' 'group_valid_RT_OE' 'group_invalid_RT_OE'};
   output_table = vertcat(output_table,...
                   horzcat(files{i}(1,1:end-4), ...
                   num2cell([size(RT,2) ...
                   behavior.group_ACC(1,i) ...
                   behavior.group_valid_ACC(1,i) ...
                   behavior.group_invalid_ACC(1,i) ...
                   behavior.group_RT(1,i) ...
                   behavior.group_valid_RT(1,i) ...
                   behavior.group_invalid_RT(1,i) ...
                   behavior.group_RT_OE(1,i) ...
                   behavior.group_valid_RT_OE(1,i) ...
                   behavior.group_invalid_RT_OE(1,i)])));
               
   writetable(cell2table(output_table, 'VariableNames', headers), fullfile(PathName, FileName));




    %Writing new EVT file
    evt_info = [time,ones(size(time,1),1),triggers];
    filename = strcat(files{i}(1,1:end-4),'_recoded.evt');
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\n',evt_header);
    fclose(fid);
    dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
end

fprintf('\nDone!\n\n')

end