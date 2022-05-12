function [tMap,total_time] = sourcespace_ancova_unpaired_DICS


%FileName #s and what they are loading in:
%1.  coh cond 1 group 1
%2.  coh cond 2 group 1
%3.  coh cond 1 group 2
%4.  coh cond 2 group 2
%5.  coh cond 1 group 3
%6.  coh cond 2 group 3
%7.  coh cond 1 group 4
%8.  coh cond 2 group 4
%
%9.  amp cond 1 group 1
%10. amp cond 2 group 1
%11. amp cond 1 group 2
%12. amp cond 2 group 2
%13. amp cond 1 group 3
%14. amp cond 2 group 3
%15. amp cond 1 group 4
%16. amp cond 2 group 4
%
%


uigetfile_message{1} = 'Select Coherence NIIs for Group 1 Condition 1';
uigetfile_message{2} = 'Select Coherence NIIs for Group 1 Condition 2';
uigetfile_message{3} = 'Select Coherence NIIs for Group 2 Condition 1';
uigetfile_message{4} = 'Select Coherence NIIs for Group 2 Condition 2';
uigetfile_message{5} = 'Select Coherence NIIs for Group 3 Condition 1';
uigetfile_message{6} = 'Select Coherence NIIs for Group 3 Condition 2';
uigetfile_message{7} = 'Select Coherence NIIs for Group 4 Condition 1';
uigetfile_message{8} = 'Select Coherence NIIs for Group 4 Condition 2';

uigetfile_message{9}  = 'Select Amplitude NIIs for Group 1 Condition 1';
uigetfile_message{10} = 'Select Amplitude NIIs for Group 1 Condition 2';
uigetfile_message{11} = 'Select Amplitude NIIs for Group 2 Condition 1';
uigetfile_message{12} = 'Select Amplitude NIIs for Group 2 Condition 2';
uigetfile_message{13} = 'Select Amplitude NIIs for Group 3 Condition 1';
uigetfile_message{14} = 'Select Amplitude NIIs for Group 3 Condition 2';
uigetfile_message{15} = 'Select Amplitude NIIs for Group 4 Condition 1';
uigetfile_message{16} = 'Select Amplitude NIIs for Group 4 Condition 2';

n_messages = length(uigetfile_message);
n_groups = 4;

for i = 1:n_messages %this is specific to 4 group DICS
    
    [FileName_temp,PathName_temp,~] = uigetfile('*.nii',uigetfile_message{i},'MultiSelect','on');
    
    FileName_list{i} = FileName_temp;
    PathName_list{i} = PathName_temp;
    
    
    cd(PathName_temp);
    % Open first NII and extract parameters%
    NII_param = load_nii(FileName_list{i}{1}); %Double {} because its a cell array of cell arrays
    if iscell(NII_param)
        NII_param = NII_param{1,1};
    end
    ref_size(i,:) = size(NII_param.img);
    nii_info = niftiinfo(FileName_list{i}{1});
    
    
    clear FileName_temp PathName_temp
    
end %end of loading for loop




if ~isequal(ref_size(1,:),ref_size(2,:),ref_size(3,:),ref_size(4,:),ref_size(5,:),ref_size(6,:),ref_size(7,:),ref_size(8,:),ref_size(9,:),ref_size(10,:),ref_size(11,:),ref_size(12,:),ref_size(13,:),ref_size(14,:),ref_size(15,:),ref_size(16,:))
    error('Volumetric images must be the same dimensions!')
elseif ~isequal(length(FileName_list{1}),length(FileName_list{9})) || ~isequal(length(FileName_list{2}),length(FileName_list{10})) || ~isequal(length(FileName_list{3}),length(FileName_list{11})) || ~isequal(length(FileName_list{4}),length(FileName_list{12})) | ~isequal(length(FileName_list{5}),length(FileName_list{13})) || ~isequal(length(FileName_list{6}),length(FileName_list{14})) || ~isequal(length(FileName_list{7}),length(FileName_list{15}))  | ~isequal(length(FileName_list{8}),length(FileName_list{16}))
    error('The same number of files must be loaded!')
end


%Collect coordinates in Talairach space and convert to double%
inputs = {'X', 'Y', 'Z'};
defaults = {'0', '0', '0'};
answer = inputdlg(inputs, 'Please Input Talairach Coordinates', 2, defaults,'on');
[X,Y,Z] = deal(answer{:});
coord_label = sprintf('X: %s; Y: %s; Z: %s',X,Y,Z);
X = str2num(X);
Y = str2num(Y);
Z = str2num(Z);

voxel_coordinates = [X,Y,Z];

%Open NII and read coordinate space parameters%
Xstart = NII_param.hdr.hist.srow_x(4);
Ystart = NII_param.hdr.hist.srow_y(4);
Zstart = NII_param.hdr.hist.srow_z(4);

voxel_coordinates(1) = voxel_coordinates(1) - Xstart;
voxel_coordinates(2) = voxel_coordinates(2) - Ystart;
voxel_coordinates(3) = voxel_coordinates(3) - Zstart;

%Convert voxel coordinates into talairach space (based on voxel size)
voxel_coordinates(1) = round(voxel_coordinates(1)/NII_param.hdr.dime.pixdim(2))+1;
voxel_coordinates(2) = round(voxel_coordinates(2)/NII_param.hdr.dime.pixdim(3))+1;
voxel_coordinates(3) = round(voxel_coordinates(3)/NII_param.hdr.dime.pixdim(4))+1;

%Check that the coordinates are within the NII space%
if voxel_coordinates(1,1) < 0 | voxel_coordinates(1,1) > size(NII_param.img,1) | voxel_coordinates(1,2) < 0 | voxel_coordinates(1,2) > size(NII_param.img,2) | voxel_coordinates(1,3) < 0 | voxel_coordinates(1,3) > size(NII_param.img,3)
    error('These seed coordinates don''t make sense! Texas Steve!');
    fprintf('\n');
end


num_per_group(1) = length(FileName_list{1});
num_per_group(2) = length(FileName_list{3});
num_per_group(3) = length(FileName_list{5});
num_per_group(4) = length(FileName_list{7});


%now you need to build the group variable column. From my testing, I think
%it NEEDS to be a cell array of strings, not numbers...
group_names = {'a','b','c','d'};

position_counter = 1;

for i = 1:length(group_names)
    
    for ii = 1:num_per_group(i)
        
        group_names_vector(position_counter,1) = group_names{i}; %not sure how to preallocate a cell array...
        
        position_counter = position_counter + 1;
        
    end
end






%Preallocate
COH_data_cond1 = zeros([(sum(num_per_group)),ref_size(1,:)]);
COH_data_cond2 = zeros([(sum(num_per_group)),ref_size(1,:)]);

AMP_data_cond1 = zeros([(sum(num_per_group)),ref_size(1,:)]);
AMP_data_cond2 = zeros([(sum(num_per_group)),ref_size(1,:)]);
AMP_data_sub   = zeros([(sum(num_per_group)),ref_size(1,:)]);


seedamp_covar_cond1 = zeros(sum(num_per_group),1);
seedamp_covar_cond2 = zeros(sum(num_per_group),1);
seedamp_covar_sub = zeros(sum(num_per_group),1);



cond1_subject_counter = 1;
cond2_subject_counter = 1;

%Load in all COH data
for i = 1:(n_groups*2)
    
    loading_counter = 1;
    
    cd(PathName_list{i})
    n_files_to_load = length(FileName_list{i});
    
    if rem(i, 2) ~= 0 %if i is odd, you are loading in condition 1
        
        while loading_counter < n_files_to_load + 1;
            
            %load in NIIs
            COH_NII = load_nii(FileName_list{i}{loading_counter});
            COH_data_cond1(cond1_subject_counter,:,:,:) = COH_NII.img;
            clear COH_NII
            
            loading_counter = loading_counter + 1;
            cond1_subject_counter = cond1_subject_counter + 1;
            
        end
    else %else, you are loading in condition 2
        
        while loading_counter < n_files_to_load + 1;
            
            %load in NIIs
            COH_NII = load_nii(FileName_list{i}{loading_counter});
            COH_data_cond2(cond2_subject_counter,:,:,:) = COH_NII.img;
            clear COH_NII
            
            loading_counter = loading_counter + 1;
            cond2_subject_counter = cond2_subject_counter + 1;
            
        end
    end
end


%Now load in amp NIIs
cond1_subject_counter = 1;
cond2_subject_counter = 1;

for i = ((n_groups*2)+1):((n_groups*2)+(n_groups*2))
    
    loading_counter = 1;
    
    cd(PathName_list{i})
    n_files_to_load = length(FileName_list{i});
    
    if rem(i, 2) ~= 0 %if i is odd, you are loading in condition 1
        
        while loading_counter < n_files_to_load + 1;
            
            %load in NIIs
            AMP_NII = load_nii(FileName_list{i}{loading_counter});
            AMP_data_cond1(cond1_subject_counter,:,:,:) = AMP_NII.img;
            clear AMP_NII
            
            seedamp_covar_cond1(cond1_subject_counter,1) = AMP_data_cond1(cond1_subject_counter,voxel_coordinates(1),voxel_coordinates(2),voxel_coordinates(3));
            
            loading_counter = loading_counter + 1;
            cond1_subject_counter = cond1_subject_counter + 1;
            
        end
    else %else, you are loading in condition 2
        
        while loading_counter < n_files_to_load + 1;
            
            %load in NIIs
            AMP_NII = load_nii(FileName_list{i}{loading_counter});
            AMP_data_cond2(cond2_subject_counter,:,:,:) = AMP_NII.img;
            clear AMP_NII
            
            seedamp_covar_cond2(cond2_subject_counter,1) = AMP_data_cond2(cond2_subject_counter,voxel_coordinates(1),voxel_coordinates(2),voxel_coordinates(3));

            
            loading_counter = loading_counter + 1;
            cond2_subject_counter = cond2_subject_counter + 1;
            
        end
    end
end


%Now subtract across the amp seed and source
AMP_data_sub = AMP_data_cond1 - AMP_data_cond2;
seedamp_covar_sub = seedamp_covar_cond1 - seedamp_covar_cond2;


clear i ii iii
F_Map_condition = zeros(ref_size(1,:));
F_Map_interaction = zeros(ref_size(1,:));
F_Map_group = zeros(ref_size(1,:));



stat_counter = 1; %Use this to print out the df text file

%Run the model
progress_bar = waitbar(0,'Performing Statistics... Please Wait...');
tic
for i = 1:ref_size(1,1)  %Looping through the x dimension
    %fprintf('Running X:%d out of %d\n',i,ref_size1(1,1))
    
    %Progress Bar
    waitbar(i/(ref_size(1,1)))
    
    for ii = 1:ref_size(1,2) %Looping through the y dimension
        for iii = 1:ref_size(1,3) %Looping through the z dimension
            if sum(COH_data_cond1(:,i,ii,iii)) == 0 || sum(COH_data_cond2(:,i,ii,iii)) == 0 || sum(AMP_data_sub(:,i,ii,iii)) == 0 
                F_Map_condition(i,ii,iii) = 0;
                F_Map_interaction(i,ii,iii) = 0;
                F_Map_group(i,ii,iii) = 0;
            elseif all([i,ii,iii] == voxel_coordinates)
                F_Map_condition(i,ii,iii) = 0;
                F_Map_interaction(i,ii,iii) = 0;
                F_Map_group(i,ii,iii) = 0;
            else


                %Make sure these variables match 
                rm_table = table(COH_data_cond1(:,i,ii,iii),...
                                 COH_data_cond2(:,i,ii,iii),...
                                 AMP_data_sub(:,i,ii,iii),...
                                 seedamp_covar_sub,...
                                 group_names_vector,...
                                 'VariableNames',...
                                 {'Cond1','Cond2','SourceAmp','SeedAmp','groups'});
                             
            
                rm_design = table([1 2]','VariableNames',{'Condition'});
                
                %repeated measures model
                rm_model = ranova(fitrm(rm_table, 'Cond1-Cond2~groups+SourceAmp+SeedAmp', 'WithinDesign', rm_design));
                
                F_Map_condition(i,ii,iii) = rm_model.F(1);
                F_Map_interaction(i,ii,iii) = rm_model.F(4);
                
                
                
                %between subjects portion
                model = anova(fitrm(rm_table, 'Cond1-Cond2~groups+SourceAmp+SeedAmp', 'WithinDesign', rm_design));

                F_Map_group(i,ii,iii) = model.F(4);



                if stat_counter == 1; %only pull df once
                    
                    df1_cond = rm_model.DF(1);
                    df2_cond = rm_model.DF(5);
                    
                    df1_group = model.DF(4);
                    df2_group = model.DF(5);
                    
                    df1_interaction = rm_model.DF(4);
                    df2_interaction = rm_model.DF(5);
                    
                    
                    stat_counter = stat_counter + 1; %Use this to print out the df text file
                    
                end

            end
        end
    end
end

cd(save_path);

df1_cond = num2str(df1_cond);
df2_cond = num2str(df2_cond);

df1_group = num2str(df1_group);
df2_group = num2str(df2_group);

df1_interaction = num2str(df1_interaction);
df2_interaction = num2str(df2_interaction);


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


