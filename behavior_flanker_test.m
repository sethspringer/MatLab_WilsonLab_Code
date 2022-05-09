function group_RT = behavior_flanker

%PURPOSE:           Calculate reaction times and statistics for the Flanker paradigm.
%
%REQUIRED INPUTS:   Cutoff Value for sd/MAD exclusions
%                       (recommend 2.5 or 3)
%                   List of evt files for calculations
%		    
%
%AUTHOR:            Brandon J. Lew, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   03/13/2019  v1: Work in Progress
%                   07/21/2021  v2: Edited by Seth Springer to output behavioral data into an excel and have a popup for outlier exclusion



%TOFIX: Some sort of accuracy exclusion
%           
%       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-Recoded Trigger Legend-%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 111 = Correct Congruent Left Stimulus 12289
% 112 = Correct Congruent Left Response 12545
% 121 = Correct Congruent Right Stimulus 12290
% 122 = Correct Congruent Right Response 12802
% 211 = Correct Incongruent Left Stimulus 12291
% 212 = Correct Incongruent Left Response 12547
% 221 = Correct Incongruent Right Stimulus 12292
% 222 = Correct Incongruent Right Response 12804
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 12289 - 12291 = stimuli
% 128* and 125* are then responses
% congruent stimuli are 12289, 90
% incongruent stimuli are 12291, 92


%%%%%%%%%%%%%%%Beginning Display Info for user%%%%%%%%%%%%%%%%%%%

DisplayText = sprintf(strcat('Welcome to the Flanker recode and behavior script!\n\n', ...
    'This script will:\n1) Read all evt files selected\n', ...
    '2) Calculate the reaction times for congruent and incongruent correct trials\n', ...
    '3) Apply stdv and mad outlier exclusions per participant\n', ...
    '4) Calculate the accuracy and average reaction time per participant\n', ...
    '5) Recode triggers and save a new evt file with the recoded triggers (see header for new codes)\n', ...
    '\nTo run this script, please and select have all of your evt files\n', ... 
    'The output of this script will be saved to your matlab workspace in a data structure'));

waitfor(msgbox(sprintf("%s\n",DisplayText)))

inputs = {'SD Cutoff'};
defaults = {'2.5'};
answer = inputdlg(inputs, 'What SD Cutoff should be used for RT outlier exclusions', 2, defaults, 'on');
[ExclusionCutoff] = deal(answer{:});
ExclusionCutoff = str2num(ExclusionCutoff);

%%Get EVTs%%

[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');  %select evt files%

evt_header = 'Tmu         	Code	TriNo';                          %create header for evt files%


cd(path);                                      %change directory to desired path

if isa(files,'char')
    iterations = 1;                        %check if only one file was selected - set iterations
else
    iterations = size(files,2);
end


%%Preallocate%%

group_RT.rawtrialdata.Cong = nan(iterations,120);    %preallocate for individual RT lists
group_RT.rawtrialdata.Incong = nan(iterations,120);
group_RT.rawtrialdata.SDexclusions.Cong = nan(iterations,50);
group_RT.rawtrialdata.SDexclusions.Incong = nan(iterations,50);
group_RT.rawtrialdata.MADexclusions.Cong = nan(iterations,50);
group_RT.rawtrialdata.MADexclusions.Incong = nan(iterations,50);

         group_RT.participant_aves.noexcl(1,1) = {'filename'};
         group_RT.participant_aves.noexcl(2,1) = {'AveCongRT'};  %saves RTs to a group matrix
         group_RT.participant_aves.noexcl(3,1) = {'AveIncongRT'};
         group_RT.participant_aves.noexcl(4,1) = {'FlankerEffect'};

         group_RT.participant_aves.MAD(1,1) = {'filename'};
         group_RT.participant_aves.MAD(2,1) = {'MADAveCongRT'};  %saves MAD excluded RTs to a group matrix
         group_RT.participant_aves.MAD(3,1) = {'MADAveIncongRT'};
         group_RT.participant_aves.MAD(4,1) = {'MADFlankerEffect'};
         
         group_RT.participant_aves.StDv(1,1) = {'filename'};
         group_RT.participant_aves.StDv(2,1) = {'sdAveCongRT'};  %saves sd excluded RTs to a group matrix
         group_RT.participant_aves.StDv(3,1) = {'sdAveIncongRT'};
         group_RT.participant_aves.StDv(4,1) = {'sdFlankerEffect'};

         group_RT.participant_aves.CoefVar(1,1) = {'filename'};
         group_RT.participant_aves.CoefVar(2,1) = {'CoefVar_Cong'};
         group_RT.participant_aves.CoefVar(3,1) = {'CoefVar_Incong'}; 

         
         
    %%%%%%%%%%%%%%%%%%% Creating .xls with behavioral data %%%%%%%%%%%%%%%%%%%%%
DefaultName = 'Flanker_Recode_Behavior_Data';
[FileName,PathName,~] = uiputfile('*.xls','Please select path for output behavioral file',DefaultName);
cd(PathName);
outputfile = fopen(FileName,'w+'); 						%open a file for writing participant data and write headers%
fprintf(outputfile, 'ParID\t group_congruent_ACC\t group_incongruent_ACC\t group_congruent_RT\t group_incongruent_RT\t group_congruent_RT_OE\t group_incongruent_RT_OE\n');        %Headers applied to the output file

         
         
         
         
         
%%Start evt for loop%%
%%
for i = 1:iterations
        if isa(files,'char')                                                                   
            evt = readBESAevt(files);          %check if only one file was selected - read evt data%
            fid = fopen(files,'rt');
        elseif length(files) > 1
            evt = readBESAevt(files{i});
            fid = fopen(files{i},'rt');
        else
            error('No compatible files selected! Sad!')
        end

        header = fgets(fid);
        fclose(fid);
       
        subjID = files{i};   %gets filenames
        group_RT.subjIDs(i,1) = {subjID};   %save subject IDs
        
        
        
        
        
        
        
    %%Start data for loop & calculate RTs%%

        numrows = size(evt); %gets size of evt
        congcounter=1;
        incongcounter=1;
        
        CongRT = 0;
        IncongRT = 0;

      for ii = 1:numrows(1)-1

         if (evt(ii,3)==12289 & (evt(ii+1,3)>12499 && evt(ii+1,3)<12600))  %if row is a button press and previous row is a stimulus
             CongRT(congcounter)=(evt(ii+1,1)-(evt(ii,1)+137000))./1000;    %compute reaction time and add to variable RT
             congcounter=congcounter+1;
             evt(ii,3)=111; %Cong, left, stimulus
             evt(ii+1,3)=112; %Cong, left, response
             
         elseif (evt(ii,3)==12290 & (evt(ii+1,3)>12799 && evt(ii+1,3)<21900))
             CongRT(congcounter)=(evt(ii+1,1)-(evt(ii,1)+137000))./1000;    %compute reaction time and add to variable RT
             congcounter=congcounter+1;
             evt(ii,3)=121; %Cong, right, stimulus
             evt(ii+1,3)=122; %Cong, right, response
             
         elseif (evt(ii,3)==12291 & (evt(ii+1,3)>12499 && evt(ii+1,3)<12600))
             IncongRT(incongcounter)=(evt(ii+1,1)-(evt(ii,1)+137000))./1000; %compute reaction time and add to variable RT
             incongcounter=incongcounter+1;
             evt(ii,3)=211; %Incong, left, stimulus
             evt(ii+1,3)=212; %Incong, left, response
             
         elseif (evt(ii,3)==12292 & (evt(ii+1,3)>12799 && evt(ii+1,3)<21900))
             IncongRT(incongcounter)=(evt(ii+1,1)-(evt(ii,1)+137000))./1000; %compute reaction time and add to variable RT
             incongcounter=incongcounter+1;
             evt(ii,3)=221; %Incong, right, stimulus
             evt(ii+1,3)=222; %Ingong, right, response

         end          
          
      end
      
          group_RT.accuracy(i,1)= congcounter-1;       %save congruent accuracy (number of congruent responses / 100 * 100) (expecting 100 trials) (times 100 for percent)
          group_RT.accuracy(i,2)= incongcounter-1;     %save incongruent accuracy (number of congruent responses / 100 * 100) (expecting 100 trials) (times 100 for percent)
      
          group_RT.rawtrialdata.Cong(i,1:length(CongRT)) = CongRT;       %save trialwise reaction times
          group_RT.rawtrialdata.Incong(i,1:length(IncongRT)) = IncongRT;   
                   
    %%End data for loop%%

    %%Reaction Time Averages and Outlier Exclusions%% 

      AveCongRT = mean(CongRT);
      CoefVar_Cong=std(CongRT)/mean(CongRT);
      CoefVar_Incong=std(IncongRT)/mean(IncongRT);

        %%Congruent MAD excludes%%
      
        lowermad = median(CongRT)-ExclusionCutoff*1.4826*mad(CongRT,1);  %calculates the lower bound for Absolute Deviation Around the Median exclusion
        uppermad = median(CongRT)+ExclusionCutoff*1.4826*mad(CongRT,1);  %calculates the upper bound for Absolute Deviation Around the Median exclusion

        CongExcludesMAD = CongRT(CongRT>uppermad | CongRT<lowermad);  %indexes all reaction times above or below the MAD exclusion bounds

        MADCongRT = CongRT;  
        MADCongRT(CongRT>uppermad | CongRT<lowermad) = [];  %Deletes all reaction times above or below the MAD exclusion bounds

        MADAveCongRT = mean(MADCongRT);  %saves new mean with MAD exclusions
        
        %%%Congruent sd excludes%%
        
        lowersd = mean(CongRT)-ExclusionCutoff*std(CongRT);  %calculates the lower bound for sdv Ave exclusion
        uppersd = mean(CongRT)+ExclusionCutoff*std(CongRT);  %calculates the upper bound for sdv Ave exclusion

        CongExcludessd = CongRT(CongRT>uppersd | CongRT<lowersd);  %indexes all reaction times above or below the sdv exclusion bounds

        sdCongRT = CongRT;  
        sdCongRT(CongRT>uppersd | CongRT<lowersd) = [];  %Deletes all reaction times above or below the sdv exclusion bounds

        sdAveCongRT = mean(sdCongRT);  %saves new mean with sdv exclusions

        AveIncongRT = mean(IncongRT);


        
        %%%Incongruent sd excludes%%
        
        lowersd = mean(IncongRT)-ExclusionCutoff*std(IncongRT);  %calculates the lower bound for sdv Ave exclusion
        uppersd = mean(IncongRT)+ExclusionCutoff*std(IncongRT);  %calculates the upper bound for sdv Ave exclusion

        IncongExcludessd = IncongRT(IncongRT>uppersd | IncongRT<lowersd);  %indexes all reaction times above or below the sdv exclusion bounds

        sdIncongRT = IncongRT;  
        sdIncongRT(IncongRT>uppersd | IncongRT<lowersd) = [];  %Deletes all reaction times above or below the sdv exclusion bounds

        sdAveIncongRT = mean(sdIncongRT);  %saves new mean with sdv exclusions



                  
        
   fprintf(outputfile, '%s\t %8.4f\t %8.4f\t %8.4f\t %8.4f\t %8.4f\t %8.4f\n',...
   files{i}(1,1:end-4),...
   group_RT.accuracy(i,1),group_RT.accuracy(i,2),...
   AveCongRT,AveIncongRT,...
   sdAveCongRT,sdAveIncongRT);        %Fields applied to the output file                



         
         
     %%Write evt%%
     
        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_recoded.evt');
        else
            filename = strcat(files{i}(1,1:end-4),'_recoded.evt');                              %save evt file%
        end
        fid = fopen(filename,'wt');
        fprintf(fid,'%s\n',evt_header);
        fclose(fid);
        dlmwrite(filename,evt,'delimiter','\t','-append','precision','%.0f');
         
       
        clear subjID evt numrows congcounter incongcounter CongRT IncongRT AveCongRT AveIncongRT FlankerEffect header 
end

%%End evt for loop%%
%%

         
%%%%%%%END MAIN FUNCTION%%%%%%%

save('behavior_flanker_result','group_RT')

%%%%%%%%%%%%%%%%%%Get current date and time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  t = now;
  c = datevec ( t );
  datetime = datestr ( c, 0 );
  
%%%%%%%%%%%%%%%%%%Write to Log File%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fclose(outputfile);

fprintf('\nDone!\n\n')

end

