function seth_sourcespace_DICS_multipleRegression

%PURPOSE:
%
%REQUIRED INPUTS:
%
%
%
%NOTES:
%
%
%AUTHOR:            Seth D. Springer, DICoN Lab, Boys Town National Research Hospital
%VERSION HISTORY:   12/5/2022  v1: First working version of program


%load the table with all subjects ParID and age...
data_table = readtable('D:\MIND_VWM\ParIDs_And_Ages.xlsx');


[FileName1,PathName1,~] = uigetfile('*.nii','Select NIIs','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});



[save_name,save_path] = uiputfile('.nii','Select the directory to save the statistical map and give a name');



n_nii = size(FileName1,2);


%Preallocate nii and seed amp arrayes
COH_data_cond1 = zeros([n_nii,ref_size1]);




%Load in NII data
cd(PathName1);
for i = 1:n_nii
    COH_NII = load_nii(FileName1{1,i});
    COH_data_cond1(i,:,:,:) = COH_NII.img;
    clear COH_NII
end


clear i


%Now grab the ages that you need by matching them to the parIDs
table_cells = table2cell(data_table);
ParIDs = table_cells(:,1);
ages   = cell2mat(table_cells(:,2));
group  = cell2mat(table_cells(:,3));

counter = 1;

for i = 1:n_nii
    
    search_sequence = FileName1{i}(1:3);
    
    if strcmp(ParIDs{counter}(1:3),search_sequence) %There are matching parIDs
        ages_OE(counter,1) = ages(counter);
        group_OE(counter,1) = group(counter);
        counter = counter + 1;
        
    else
        
        while ~(strcmp(ParIDs{counter}(1:3),search_sequence))
            
            counter = counter + 1;
            
        end
        
        if strcmp(ParIDs{counter}(1:3),search_sequence) %There are matching parIDs
            ages_OE(counter,1) = ages(counter);
            group_OE(counter,1) = group(counter);
            counter = counter + 1;
            
        end
    end
end

%Remove zeros bc these are missing participants
ages_OE(ages_OE==0) = [];
group_OE(group_OE==0) = [];



clear i
F_Map_interaction = zeros(ref_size1);

stat_counter = 1; %Use this to print out the df text file
%Run the model

progress_bar = waitbar(0,'Performing Statistics... Please Wait...');

tic

for i = 1:ref_size1(1,1)  %Looping through the x dimension
    %fprintf('Running X:%d out of %d\n',i,ref_size1(1,1))
    
    %Progress Bar
    waitbar(i/(ref_size1(1,1)))
    
    for ii = 1:ref_size1(1,2) %Looping through the y dimension
        for iii = 1:ref_size1(1,3) %Looping through the z dimension
            if sum(COH_data_cond1(:,i,ii,iii)) == 0  %can't run stats on zeros
                                
            else
                
                
                
                
                
                %Creating the variables table
                variable_table = table(...
                    COH_data_cond1(:,i,ii,iii),...
                    group_OE,...
                    ages_OE,...
                    'VariableNames',...
                    {'NII','group','age'});
                
                model = fitlm(variable_table, 'NII~group*age');
                
               

                %The position of the F-value of interest depends on the order that factors were entered into the model
                F_Map_interaction(i,ii,iii) =  model.Coefficients.tStat(4);
                
                if stat_counter == 1; %only pull df once
                    
                    df1 = model.DFE;
                    
                    
                    stat_counter = stat_counter + 1; %Use this to print out the df text file
                    
                end
                
                
            end
        end
    end
end

cd(save_path);

df1 = num2str(df1);
df2 = num2str(df2);


close(progress_bar)


save_name = strcat(save_name(1:end-4),'_DF_',df1,'_',df2);

%Save the condition contrast F-map
niftiwrite(F_Map_interaction,save_name,nii_info);

total_time = toc;


fprintf('\nThese stats took %4.2f seconds in total\n\n',total_time)

end

