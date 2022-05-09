function PESCAH_OneDot_OlderPar_Recode(trials_to_remove)

%PURPOSE:           Randomly removes "trials_to_remove" trials from these subjects.
%
%REQUIRED INPUTS:   EVT files for recoding.
%
%		    
%
%NOTES:             Changes non-starting 12309s to 20
%                   Changes removed starting 12309s to 30s
%		   
%
%                  
%AUTHOR:            Seth, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   06/28/2021  v1: First working version of program

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

    
cd(path);  


if isa(files,'char')                                %set number of iterations based on whether there is one or multiple files%
    iter = 1;
else
    iter = length(files);
end


    evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%



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
    
    
    %Loop through each person's evt
    for ii = 2:length(triggers);
        if triggers(ii) == 12309;
            prev_trig = triggers(ii-1);
            if prev_trig == 12343 | prev_trig == 12351;
                triggers(ii) = 12309; %Trial start
           
            else
                triggers(ii) = 20; %12309 that isn't a trial start
            end
        end
    end
    

    
    %Randomly recode some starting triggers
    entrainment_start_index = find(triggers==12309);
    
    number_trials = length(entrainment_start_index);
    
    rand_num_vector = randperm(number_trials);
    
    subset_rand_num = rand_num_vector(1:trials_to_remove);
    
    
    positions_to_delete = entrainment_start_index(subset_rand_num);
 
    
    
    %Recode the starting triggers that you randomly picked.
    for iii = 1:length(positions_to_delete)
        triggers(positions_to_delete(iii)) = 30;
    end
    
    
    
   %Stitch an EVT back together
        evt_info = [time,data(1:size(time,1),2),triggers];
        
        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_recoded_removeTrials.evt');
        else
            filename = strcat(files{i}(1,1:end-4),'_recoded_removeTrials.evt');                              %save evt file%
        end
        
        fid = fopen(filename,'wt');
        fprintf(fid,'%s\n',evt_header);
        fclose(fid);
        dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
        
end

end