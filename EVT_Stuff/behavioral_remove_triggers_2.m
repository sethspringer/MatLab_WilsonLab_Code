function behavioral_remove_trigger(numberTrialsToRemove)

%PURPOSE:           Randomly remove a percentage of a certain trigger.
%
%REQUIRED INPUTS:   EVT files.
%
%		    
%
%NOTES:            Finds all "triggerNumber" triggers
%		           Randomly removes "percentageToRemove" of them.
%
%
%                  
%AUTHOR:            Seth D. Springer, IHN
%VERSION HISTORY:   05/26/2021  v1: First working version of program


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%
    


cd(path);  
if isa(files,'char')                                                                        %set number of iterations based on whether there is one or multiple files%
    iter = 1;
else
    iter = length(files);
end



triggerNumber1 = 31;
triggerNumber2 = 32;


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

    whereTriggerNumber = find(triggers==triggerNumber1 | triggers==triggerNumber2);


    randomTriggerNumberVector = whereTriggerNumber(randperm(length(whereTriggerNumber)));


    shortenedRandomVector = randomTriggerNumberVector(1:numberTrialsToRemove)


    triggers(shortenedRandomVector)=999

    

        evt_info = [time,data(1:size(time,1),2),triggers];
        if isa(files,'char')
            filename = strcat(files(1,1:end-4),'_removedTriggers.evt');
        else
            filename = strcat(files{i}(1,1:end-4),'_removedTriggers.evt');                              %save evt file%
        end
        fid = fopen(filename,'wt');
        fprintf(fid,'%s\n',evt_header);
        fclose(fid);
        dlmwrite(filename,evt_info,'delimiter','\t','-append','precision','%.0f');
end