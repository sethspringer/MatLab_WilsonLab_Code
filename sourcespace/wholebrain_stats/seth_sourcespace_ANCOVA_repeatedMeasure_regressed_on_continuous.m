function seth_sourcespace_ANCOVA_repeatedMeasure_regressed_on_continuous

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
%VERSION HISTORY:   4/21/2022  v1: First working version of program


%load the table with all subjects ParID and age...
ages_table = readtable('D:\MIND_VWM\ParIDs_And_Ages.xlsx');


[FileName1,PathName1,~] = uigetfile('*.nii','Select NIIs for Condition 1 Group 1','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});

[FileName2,PathName2,~] = uigetfile('*.nii','Select NIIs for Condition 2 Group 1','MultiSelect','on');
cd(PathName2);
% Open first NII and extract parameters%
NII_param2 = load_nii(FileName2{1,1});
if iscell(NII_param2)
    NII_param2 = NII_param2{1,1};
end
ref_size2 = size(NII_param2.img);



if ~isequal(ref_size1,ref_size2)
    error('Volumetric images must be the same dimensions!')
elseif ~isequal(size(FileName1),size(FileName2))
    error('The same number of files must be loaded for each condition for Group 1!')
end


save_path = uigetdir(PathName1,'Select the directory to save the statistical maps');



n_nii = size(FileName1,2);



COH_data_cond1 = zeros([n_nii,ref_size1]);
COH_data_cond2 = zeros([n_nii,ref_size1]);



%Load in condition 1 data
cd(PathName1);
for i = 1:n_nii
    COH_NII = load_nii(FileName1{1,i});
    COH_data_cond1(i,:,:,:) = COH_NII.img;
    clear COH_NII
end


clear i

%Load in condition 2 data
cd(PathName2);
for i = 1:n_nii
    COH_NII = load_nii(FileName2{1,i});
    COH_data_cond2(i,:,:,:) = COH_NII.img;
    clear COH_NII
end





%Now grab the ages that you need by matching them to the parIDs
ages_cells = table2cell(ages_table);
ParIDs = ages_cells(:,1);
ages   = cell2mat(ages_cells(:,2));

counter = 1;

for i = 1:n_nii
    
    search_sequence = FileName1{i}(1:3);
    
    if strcmp(ParIDs{counter}(1:3),search_sequence) %There are matching parIDs
        ages_OE(counter,1) = ages(counter);
        counter = counter + 1;
        
    else
        
        while ~(strcmp(ParIDs{counter}(1:3),search_sequence))
        
        counter = counter + 1;
        
        end
        
        if strcmp(ParIDs{counter}(1:3),search_sequence) %There are matching parIDs
            ages_OE(counter,1) = ages(counter);
            counter = counter + 1;
            
        end
        
    end
    
end

%Remove zeros bc these are missing participants
ages_OE(ages_OE==0) = [];


clear i
F_Map_condition = zeros(ref_size1);
F_Map_interaction = zeros(ref_size1);
F_Map_age = zeros(ref_size1);

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
            if sum(COH_data_cond1(:,i,ii,iii)) == 0 || sum(COH_data_cond2(:,i,ii,iii)) == 0

            else
                %Creating the variables table
                rm_table = table(COH_data_cond1(:,i,ii,iii),...
                                 COH_data_cond2(:,i,ii,iii),...
                                 ages_OE,...
                                 'VariableNames',...
                                 {'Cond1','Cond2','age'});
            
                rm_design = table([1 2]','VariableNames',{'Condition'});;
                rm_model = ranova(fitrm(rm_table, 'Cond1-Cond2~age', 'WithinDesign', rm_design));
                
                
                %The position of the F-value of interest depends on the order that factors were entered into the model
                F_Map_condition(i,ii,iii) = rm_model.F(1);
                F_Map_interaction(i,ii,iii) = rm_model.F(2);
                
                if stat_counter == 1; %only pull df once
                    
                    df1 = rm_model.DF(1);
                    df2 = rm_model.DF(3);
                    stat_counter = stat_counter + 1; %Use this to print out the df text file
                    
                end

                
                %Now do the main effect of age
                COH_both = [COH_data_cond1(:,i,ii,iii),COH_data_cond2(:,i,ii,iii)];
                COH_ave = mean(COH_both,2);
                
                rm_table_age_effect = table(COH_ave,...
                                              ages_OE,...
                                              'VariableNames',...
                                              {'Cond','age'});
                
                
                age_effect_model = fitlm(rm_table_age_effect,'Cond~age');
                
                F_Map_age(i,ii,iii) = (age_effect_model.Coefficients.tStat(2))^2; %Pull the t-value and sqr it to get the F-value
                
                
            end
        end
    end
end

cd(save_path);

df1 = num2str(df1);
df2 = num2str(df2);


close(progress_bar)


condition_file_name = strcat('ConditionContrast_DF_',df1,'_',df2);
age_file_name = strcat('AgeContrast_DF_',df1,'_',df2);
interaction_file_name = strcat('GroupByCondition_Interaction_DF_',df1,'_',df2);




%Save the condition contrast F-map
niftiwrite(F_Map_condition,condition_file_name,nii_info);

%Save the group contrast F-map
niftiwrite(F_Map_age,age_file_name,nii_info);

%Save the interaction F-map
niftiwrite(F_Map_interaction,interaction_file_name,nii_info);

total_time = toc;


fprintf('\nThese stats took %4.2f seconds in total\n\n',total_time)

end

