function [behavior] = VWM_behavior_general_recode

%PURPOSE:           Computes behavior for VWM.
%
%REQUIRED INPUTS:   EVT files 
%
%		    
%
%NOTES:            Removes ALL 12295 triggers; removes
%                   last 12289/12290 
%                  Recodes all responses: 128**->12800, 125**->12545
%                  Computes behavior	   
%
%
%                  
%AUTHOR:            Christine M. Embury, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   04/20/20  v1: First working version of program;
%                                 modified from posner_recode_v2 Alex I Wiesman code
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixation, encoding, maintenance, retrieval, button presses

%12292 fixation
%12291 Encoding
%12293 Maintenance
%12295 nuisance following 93 always
%12289/12290 retrieval
%12289->125*; 12290->128*
%Recode correct hit to 31, correct negative to 32
%Recoded incorrect trials as well: Miss to 61, False Alarm to 62 


inputs = {'SD Cutoff'};
defaults = {'2.5'};
answer = inputdlg(inputs, 'What SD Cutoff should be used for w/in-subject RT outlier exclusions', 2, defaults, 'on');
[sd_cutoff] = deal(answer{:});
sd_cutoff = str2num(sd_cutoff);


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%
    
cd(path);  
if isa(files,'char')                                                                        %set number of iterations based on whether there is one or multiple files%
    iter = 1;
else
    iter = length(files);
end

behavior.group_ACC = NaN(1,iter);

behavior.group_RT = NaN(1,iter);
behavior.group_correct_RT = NaN(1,iter);
behavior.group_incorrect_RT = NaN(1,iter);

behavior.group_RT_OE = NaN(1,iter);
behavior.group_correct_RT_OE = NaN(1,iter);
behavior.group_incorrect_RT_OE = NaN(1,iter);

behavior.group_correct_trial_count = NaN(1,iter);
behavior.group_incorrect_trial_count = NaN(1,iter);


behavior.group_dprime = NaN(1,iter);




for i = 1:iter                                                                              %loop for files%
    if isa(files,'char')                                                                    %read evt data for file i%
        data = readBESAevt(files);
    elseif length(files) > 1;
        data = readBESAevt(files{i});
    else
        error('No compatible files selected! Boom roasted!')
    end
    
    triggers = data(:,3);
    time = data(:,1);
    time(triggers==12295) = [];                     %delete nuisance trigger 12295 from data%
    triggers(triggers==12295) = [];                     
    time(triggers==12292) = [];                     %delete fixation trigger 12292 from data%
    triggers(triggers==12292) = []; 
    
    
    
    triggers((12800<triggers)&(triggers<12810)) = 12800;    %recode the '12800-12810' triggers as '12800'
    triggers((12540<triggers)&(triggers<12555)) = 12540;    %recode the '12540-12555' triggers as '12540'
    if triggers(length(triggers))==12289|12290
        triggers(length(triggers))=12288;
    end
    

    
    
    
    
    
%     for ii = 3:length(triggers)
%         if (triggers(ii) == 12289) && (triggers(ii+1) == 12540)% HIT; correct trial%
%                 triggers(ii)=31;
%         elseif (triggers(ii)==12290)&& (triggers(ii+1)== 12800) %CORRECT NEGATIVE; correct trial%
%                 triggers(ii)=32;
%         elseif (triggers(ii) == 12289) && (triggers(ii+1) == 12800) %MISS; incorrect trial%
%                 triggers(ii)=61;
%         elseif (triggers(ii)==12290)&& (triggers(ii+1)== 12540) %FALSE ALARM; incorrect trial%
%                 triggers(ii)=62;
%         end
%     end
%     
    
      for ii = 3:length(triggers)
         if  (triggers(ii) == 12289) && (triggers(ii+1) == 12540)% HIT; correct trial%
                triggers(ii)=31;
         end
      end   
    
    
    for ii = 3:length(triggers)
         if (triggers(ii)==12290)&& (triggers(ii+1)== 12800) %CORRECT NEGATIVE; correct trial%
                triggers(ii)=32;
         end
     end
     
     for ii = 3:length(triggers)
         if (triggers(ii) == 12289) && (triggers(ii+1) == 12800) %MISS; incorrect trial%
                triggers(ii)=61;
         end
     end  
    
     for ii = 3:length(triggers)
         if (triggers(ii)==12290)&& (triggers(ii+1)== 12540) %FALSE ALARM; incorrect trial%
                triggers(ii)=62;
         end
     end  
     
    
     
    time(triggers==0) = [];                     
    triggers(triggers==0) = [];
    
    
%     for ii = 3:length(triggers)
%         if triggers(ii) == 12291
%             prev_trig = num2str(triggers(ii-1));
%             if strcmp(prev_trig(1:3),'125') | strcmp(prev_trig(1:3),'128') | strcmp(prev_trig,'12289') | strcmp(prev_trig,'12290')%                                                     %keep these%
%             else
%                 triggers(ii) = NaN;
%                 time(iial_counter = 0;
%
%             end
%         end
%     end
    triggers(isnan(triggers)) = [];
    time(isnan(time)) = [];
    
%     correct_trial_counter = 0;
%     incorrect_trial_counter = 0;
 
    RT = NaN(1,length(triggers));
    trial_index = NaN(1,length(triggers));
    trial_index_dprime = NaN(1,length(triggers));
    for ii = 3:length(triggers)
        if triggers(ii) == 31% HIT; correct trial%
                trial_index_dprime(1,ii)= 1;
                trial_index(1,ii) = 1;
                RT(1,ii) = ((time(ii+1)-time(ii))-33000)/1000;
        elseif triggers(ii)==32 %CORRECT NEGATIVE; correct trial%
                trial_index_dprime(1,ii)= 2;
                trial_index(1,ii) = 1;
                RT(1,ii) = ((time(ii+1)-time(ii))-33000)/1000;            
        elseif triggers(ii) == 61 %MISS; incorrect trial%
             trial_index_dprime(1,ii)= 3;
             trial_index(1,ii) = 0;
             RT(1,ii) = ((time(ii+1)-time(ii))-33000)/1000;
        elseif triggers(ii)==62 %FALSE ALARM; incorrect trial%
             trial_index_dprime(1,ii)= 4;
             trial_index(1,ii) = 0;
             RT(1,ii) = ((time(ii+1)-time(ii))-33000)/1000;
        end
    end
    
    RT(isnan(RT)) = [];
    trial_index(isnan(trial_index)) = []; 
    trial_index_dprime(isnan(trial_index_dprime))=[];
    
%     if isa(files,'char')
%         fprintf('Recoded %d trials for file %s\n',size(RT,2),files);
%     else
%         fprintf('Recoded %d trials for file %s\n',size(RT,2),files{i});
%     end
%     
    if size(RT,2) > 1 && size(trial_index,2) > 1
        %behavior.group_ID(1,i) = 
        behavior.group_ACC(1,i) = ((sum(trial_index))/128)*100;
        %behavior.group_ACC(1,i) = size(validity_index,2)/(correct_trial_counter+incorrect_trial_counter);
        %behavior.group_correct_ACC(1,i) = sum(validity_index == 1)/correct_trial_counter;
        %behavior.group_incorrect_ACC(1,i) = sum(validity_index == 0)/incorrect_trial_counter;
        
        behavior.group_correct_trial_count(1,i)=sum(trial_index);
        behavior.group_incorrect_trial_count(1,i)=sum(~trial_index);
        
        trial_index = logical(trial_index);
        behavior.group_RT(1,i) = mean(RT,2);
        behavior.group_correct_RT(1,i) = mean(RT(trial_index),2);
        behavior.group_incorrect_RT(1,i) = mean(RT(~trial_index),2);
        
        RT_OE = RT(RT<(mean(RT)+(sd_cutoff*std(RT))) & RT>(mean(RT)-(sd_cutoff*std(RT))));
        trial_index_OE = trial_index(RT<(mean(RT)+(sd_cutoff*std(RT))) & RT>(mean(RT)-(sd_cutoff*std(RT))));
        behavior.group_RT_OE(1,i) = mean(RT_OE,2);
        behavior.group_correct_RT_OE(1,i) = mean(RT_OE(trial_index_OE),2);
        behavior.group_incorrect_RT_OE(1,i) = mean(RT_OE(~trial_index_OE),2);
    
       

        
       behavior.group_dprime(1,i) = norminv(((sum(trial_index_dprime==1))/(sum(trial_index_dprime==1)+sum(trial_index_dprime==3))))-norminv(((sum(trial_index_dprime==4))/(sum(trial_index_dprime==4)+sum(trial_index_dprime==2))));
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