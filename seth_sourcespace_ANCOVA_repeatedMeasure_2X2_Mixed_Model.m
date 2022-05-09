function seth_sourcespace_ANCOVA_repeatedMeasure_2X2_Mixed_Model

%PURPOSE:           Perform a 2x2 mixed model (group-by-condition) contrast on DICS
%                   connectivity maps, controlling for amplitude at the source voxel and seed.
%
%REQUIRED INPUTS:   DICS NIIs for each group and condition, amplitude maps that are subtracted condition-wise, 
%                   and the seed DICS voxel location in TAL space coordinates
%
%		    
%
%NOTES:            
%		  
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, Boys Town National Research Hospital
%VERSION HISTORY:   3/10/2022  v1: First working version of program



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

[FileName3,PathName3,~] = uigetfile('*.nii','Select NIIs for Condition 1 Group 2','MultiSelect','on');
cd(PathName3);
% Open first NII and extract parameters%
NII_param3 = load_nii(FileName3{1,1});
if iscell(NII_param3)
    NII_param3 = NII_param3{1,1};
end
ref_size3 = size(NII_param3.img);


[FileName4,PathName4,~] = uigetfile('*.nii','Select NIIs for Condition 2 Group 2','MultiSelect','on');
cd(PathName4);
% Open first NII and extract parameters%
NII_param4 = load_nii(FileName4{1,1});
if iscell(NII_param2)
    NII_param4 = NII_param4{1,1};
end
ref_size4 = size(NII_param4.img);




if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4)
    error('Volumetric images must be the same dimensions!')
elseif ~isequal(size(FileName1),size(FileName2))
    error('The same number of files must be loaded for each condition for Group 1!')
elseif ~isequal(size(FileName3),size(FileName4))
    error('The same number of files must be loaded for each condition for Group 2!')
end


save_path = uigetdir(PathName1,'Select the directory to save the statistical maps');



n_group1 = size(FileName1,2);
n_group2 = size(FileName3,2);
n_both_groups = n_group1+n_group2;


COH_data_cond1 = zeros([n_both_groups,ref_size1]);
COH_data_cond2 = zeros([n_both_groups,ref_size1]);



%Load in condition 1 data
cd(PathName1);
cond1_counter = 1; %This counter is needed bc data is loaded in in two separate for loops
for i = 1:n_group1
    COH_NII = load_nii(FileName1{1,i});
    COH_data_cond1(cond1_counter,:,:,:) = COH_NII.img;
    cond1_counter = cond1_counter + 1;
    clear COH_NII
end

cd(PathName3);
for i = 1:n_group2
    COH_NII = load_nii(FileName3{1,i});
    COH_data_cond1(cond1_counter,:,:,:) = COH_NII.img;
    cond1_counter = cond1_counter + 1;    
    clear COH_NII
end

%Load in condition 2 data
cd(PathName2);
cond2_counter = 1; %This counter is needed bc data is loaded in in two separate for loops
for i = 1:n_group1
    COH_NII = load_nii(FileName2{1,i});
    COH_data_cond2(cond2_counter,:,:,:) = COH_NII.img;
    cond2_counter = cond2_counter + 1;
    clear COH_NII
end

cd(PathName4);
for i = 1:n_group2
    COH_NII = load_nii(FileName4{1,i});
    COH_data_cond2(cond2_counter,:,:,:) = COH_NII.img;
    cond2_counter = cond2_counter + 1;    
    clear COH_NII
end



%need to create a group variable of -1s and 1s.

group1_contrast_code_array = (ones(n_group1,1))*-1; %contrast codes for group 1 is -1
group2_contrast_code_array = ones(n_group2,1); %contrast codes for group 2 is 1
group = [group1_contrast_code_array;group2_contrast_code_array]; %one contrast code array for the model


clear i
F_Map_condition = zeros(ref_size1);
F_Map_interaction = zeros(ref_size1);
F_Map_group = zeros(ref_size1);

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
                F_Map_condition(i,ii,iii) = 0;
                F_Map_interaction(i,ii,iii) = 0;
            else
                %Creating the variables table
                rm_table = table(COH_data_cond1(:,i,ii,iii),...
                                 COH_data_cond2(:,i,ii,iii),...
                                 group,...
                                 'VariableNames',...
                                 {'Cond1','Cond2','group'});
            
                rm_design = table([1 2]','VariableNames',{'Condition'});;
                rm_model = ranova(fitrm(rm_table, 'Cond1-Cond2~group', 'WithinDesign', rm_design));
                
                
                %The position of the F-value of interest depends on the order that factors were entered into the model
                F_Map_condition(i,ii,iii) = rm_model.F(1);
                F_Map_interaction(i,ii,iii) = rm_model.F(2);
                
                if stat_counter == 1; %only pull df once
                    
                    df1 = rm_model.DF(1);
                    df2 = rm_model.DF(3);
                    stat_counter = stat_counter + 1; %Use this to print out the df text file
                    
                end

                
                %Now do the main effect of group
                COH_both = [COH_data_cond1(:,i,ii,iii),COH_data_cond2(:,i,ii,iii)];
                COH_ave = mean(COH_both,2);
                
                rm_table_group_effect = table(COH_ave,...
                                              group,...
                                              'VariableNames',...
                                              {'Cond','group'});
                
                
                group_effect_model = fitlm(rm_table_group_effect,'Cond~group');
                
                F_Map_group(i,ii,iii) = (group_effect_model.Coefficients.tStat(2))^2; %Pull the t-value and sqr it to get the F-value
                
                
            end
        end
    end
end

cd(save_path);

df1 = num2str(df1);
df2 = num2str(df2);


close(progress_bar)


condition_file_name = strcat('ConditionContrast_DF_',df1,'_',df2);
group_file_name = strcat('GroupContrast_DF_',df1,'_',df2);
interaction_file_name = strcat('GroupByCondition_Interaction_DF_',df1,'_',df2);




%Save the condition contrast F-map
niftiwrite(F_Map_condition,condition_file_name,nii_info);

%Save the group contrast F-map
niftiwrite(F_Map_group,group_file_name,nii_info);

%Save the interaction F-map
niftiwrite(F_Map_interaction,interaction_file_name,nii_info);

total_time = toc;


fprintf('\nThese stats took %4.2f seconds in total\n\n',total_time)

end

