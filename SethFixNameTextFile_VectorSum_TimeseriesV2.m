function VectorSum_Timeseries(average)

%PURPOSE:           Compute the vector sum of and then concatenate peak voxel
%                   timeseries.
%
%REQUIRED INPUTS:   average = 'y' (to also return an averaged timecourse of
%                   each orientation/orientation combination) or 'n'
%
%NOTES:             (1)  Function will read any and all ".tfc" files within a directory with
%                        "Timeseries" in their title
%                   (2)  Function will read files in numerical order according to their titles!
%                  
%AUTHOR:            Alex I. Wiesman, University of Nebraska Medical Center
%VERSION HISTORY:   04/07/2017  v1: First working version of program

%Edit V2: Gives a FileName.txt with all the filenames, Brandon Lew 1/21/2020


if average == 'What is the meaning of life?'
    fprintf('It''s all about that Mulan Szechuan McNugget Sauce Morty!!!\n')
else
%--------------------------------------------------------------%
% -- Compute Vector Sums and Concatenate All TFCs Vertically --%
%--------------------------------------------------------------%

files = dir('*Timeseries*.tfc');
number = length(files);
FileNameConcat = {};
if number == 0
    fprintf('No TFCs found! Ya done goofed! \n');
else
    for i = 1
        FileName = files(i).name;
        TFC = dlmread(FileName);
        vectorsum = sqrt((TFC(1,:).^2)+(TFC(2,:).^2));
        TFCconcat = vectorsum;
        FileNameConcat = [FileNameConcat;FileName];    %I ONLY ADDED THIS LINE!
        for k = 2:length(files)
            FileName = files(k).name;
            TFC = dlmread(FileName);
            vectorsum = sqrt((TFC(1,:).^2)+(TFC(2,:).^2));
            TFCconcat = [TFCconcat;vectorsum];
            FileNameConcat = [FileNameConcat;FileName];
        end
    end
end

%-----------------------------------------------------------------------%
%-- Order Vector Sum Timeseries and (Optionally) Average By Timepoint --%
%-----------------------------------------------------------------------%

    if average == 'n'
        dlmwrite('VectorSum_Timeseries.txt',TFCconcat,'delimiter','\t','precision',4) 
    elseif average == 'y'
        dlmwrite('VectorSum_Timeseries.txt',TFCconcat,'delimiter','\t','precision',4)
        AveragedVectorSums = mean(TFCconcat);
        dlmwrite('VectorSum_Timeseries_Means.txt',AveragedVectorSums,'delimiter','\t','precision',4)
    
    end
    
FileNameTable=cell2table(FileNameConcat);
writetable(FileNameTable,'FileNames.txt');
    
fprintf('%.0f files summed and concatenated! Savage! \n', number);
end
end