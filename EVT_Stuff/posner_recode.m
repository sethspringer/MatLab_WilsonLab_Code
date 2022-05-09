function posner_recode

%PURPOSE:           Recode triggers for posner.
%
%REQUIRED INPUTS:   EVT files for recoding MAMMOA posner triggers.
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


[files,path,~] = uigetfile('*.evt','Please Select EVTs','Multiselect','on');        %select evt files%

evt_header = 'Tmu         	Code	TriNo';                                         %create header for evt files%
    
cd(path);  
if isa(files,'char')                                                                        %set number of iterations based on whether there is one or multiple files%
    iter = 1;
else
    iter = length(files);
end

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
    time(triggers==12303) = [];                     %delete nonsense triggers 12303 and 12302 from data%
    triggers(triggers==12303) = [];
    time(triggers==12302) = [];                     
    triggers(triggers==12302) = [];
    
    time(triggers==0) = [];                     
    triggers(triggers==0) = [];
    
    for ii = 3:length(triggers)
        if triggers(ii) == 12299
            prev_trig = num2str(triggers(ii-1));
            if strcmp(prev_trig(1:3),'125') | strcmp(prev_trig(1:3),'128') | strcmp(prev_trig,'12289') | strcmp(prev_trig,'12290') | strcmp(prev_trig,'12291') | strcmp(prev_trig,'12292') | strcmp(prev_trig,'12293') | strcmp(prev_trig,'12294') | strcmp(prev_trig,'12295') | strcmp(prev_trig,'12296')
                                                    %keep these%
            else
            triggers(ii) = 11111;
            time(ii) = 11111;
            end
        elseif triggers(ii) == 12295
            prev_trig = num2str(triggers(ii-1));
            if strcmp(prev_trig,'12297') | strcmp(prev_trig,'12298')
                                                    %keep these%
            else
            triggers(ii) = 11111;
            time(ii) = 11111;
            end
        end
    end
  triggers(triggers==11111) = [];
  time(time==11111) = [];
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