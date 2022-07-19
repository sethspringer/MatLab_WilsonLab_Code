function seth_sourcespace_DICS_ANCOVA_repeatedMeasure_oneWay_3condition

%PURPOSE:           
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
%VERSION HISTORY:   7/15/2022  v1: First working version of program



[FileName1,PathName1,~] = uigetfile('*.nii','Select Coherence NIIs for Condition 1','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});

[FileName2,PathName2,~] = uigetfile('*.nii','Select Coherence NIIs for Condition 2','MultiSelect','on');
cd(PathName2);
% Open first NII and extract parameters%
NII_param2 = load_nii(FileName2{1,1});
if iscell(NII_param2)
    NII_param2 = NII_param2{1,1};
end
ref_size2 = size(NII_param2.img);

[FileName3,PathName3,~] = uigetfile('*.nii','Select Coherence NIIs for Condition 3','MultiSelect','on');
cd(PathName3);
% Open first NII and extract parameters%
NII_param3 = load_nii(FileName3{1,1});
if iscell(NII_param3)
    NII_param3 = NII_param3{1,1};
end
ref_size3 = size(NII_param3.img);


[FileName4,PathName4,~] = uigetfile('*.nii','Select Amp NIIs for Condition 1','MultiSelect','on');
cd(PathName4);
% Open first NII and extract parameters%
NII_param4 = load_nii(FileName4{1,1});
if iscell(NII_param2)
    NII_param4 = NII_param4{1,1};
end
ref_size4 = size(NII_param4.img);



[FileName5,PathName5,~] = uigetfile('*.nii','Select Amp NIIs for Condition 2','MultiSelect','on');
cd(PathName5);
% Open first NII and extract parameters%
NII_param5 = load_nii(FileName5{1,1});
if iscell(NII_param5)
    NII_param5 = NII_param5{1,1};
end
ref_size5 = size(NII_param5.img);


[FileName6,PathName6,~] = uigetfile('*.nii','Select Amp NIIs for Condition 3','MultiSelect','on');
cd(PathName6);
% Open first NII and extract parameters%
NII_param6 = load_nii(FileName6{1,1});
if iscell(NII_param6)
    NII_param6 = NII_param6{1,1};
end
ref_size6 = size(NII_param6.img);



if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4,ref_size5,ref_size6)
    error('Volumetric images must be the same dimensions!')
elseif ~isequal(size(FileName1),size(FileName2),size(FileName3),size(FileName4),size(FileName5),size(FileName6))
    error('The same number of files must be loaded for each condition and the amplitude covariates!')
end


%Collect coordinates in Talairach space and convert to double%

inputs = {'X', 'Y', 'Z'};
defaults = {'-34', '-28.5', '49.5'};	
answer = inputdlg(inputs, 'Please Input Talairach Coordinates', 2, defaults,'on');
[X,Y,Z] = deal(answer{:});
coord_label = sprintf('X: %s; Y: %s; Z: %s',X,Y,Z);
X = str2num(X);
Y = str2num(Y);
Z = str2num(Z);

voxel_coordinates = [X,Y,Z];

%Open NII and read coordinate space parameters%
Xstart = NII_param1.hdr.hist.srow_x(4);
Ystart = NII_param1.hdr.hist.srow_y(4);
Zstart = NII_param1.hdr.hist.srow_z(4);

voxel_coordinates(1) = voxel_coordinates(1) - Xstart;
voxel_coordinates(2) = voxel_coordinates(2) - Ystart;
voxel_coordinates(3) = voxel_coordinates(3) - Zstart;

%Convert voxel coordinates into talairach space (based on voxel size)
voxel_coordinates(1) = round(voxel_coordinates(1)/NII_param1.hdr.dime.pixdim(2))+1;
voxel_coordinates(2) = round(voxel_coordinates(2)/NII_param1.hdr.dime.pixdim(3))+1;
voxel_coordinates(3) = round(voxel_coordinates(3)/NII_param1.hdr.dime.pixdim(4))+1;

%Check that the coordinates are within the NII space%
if voxel_coordinates(1,1) < 0 | voxel_coordinates(1,1) > size(NII_param1.img,1) | voxel_coordinates(1,2) < 0 | voxel_coordinates(1,2) > size(NII_param1.img,2) | voxel_coordinates(1,3) < 0 | voxel_coordinates(1,3) > size(NII_param1.img,3)
    error('These seed coordinates don''t make sense! Texas Steve!');
	fprintf('\n');
end


%Where to save output files
%[save_name_cond,save_path_cond,save_index_cond] = uiputfile('*.nii','Please save condition contrast NII');
%[save_name_interaction,save_path_interaction,save_index_interaction] = uiputfile('*.nii','Please save condition-by-group interaction NII');

save_path = uigetdir(PathName1,'Select the directory to save the statistical maps');



n_group1 = length(FileName1);


COH_data_cond1 = zeros([n_group1,ref_size1]);
COH_data_cond2 = zeros([n_group1,ref_size1]);
COH_data_cond3 = zeros([n_group1,ref_size1]);

AMP_data_cond1 = zeros([n_group1,ref_size1]);
AMP_data_cond2 = zeros([n_group1,ref_size1]);
AMP_data_cond3 = zeros([n_group1,ref_size1]);



seedamp_covar = zeros(n_group1,3);


%Load in condition 1 data
cd(PathName1);
for i = 1:n_group1
    COH_NII = load_nii(FileName1{1,i});
    COH_data_cond1(i,:,:,:) = COH_NII.img;
    clear COH_NII
end

%Load in condition 2 data
cd(PathName2);
for i = 1:n_group1
    COH_NII = load_nii(FileName2{1,i});
    COH_data_cond2(i,:,:,:) = COH_NII.img;
    clear COH_NII
end

%Load in condition 3 data
cd(PathName3);
for i = 1:n_group1
    COH_NII = load_nii(FileName3{1,i});
    COH_data_cond3(i,:,:,:) = COH_NII.img;
    clear COH_NII
end



%Load in the amplitude data for condition 1
cd(PathName4);
for i = 1:n_group1
    AMP_NII = load_nii(FileName4{1,i});
    AMP_data_cond1(i,:,:,:) = AMP_NII.img;
    
    seedamp_covar(i,1) = AMP_data_cond1(i,voxel_coordinates(1), voxel_coordinates(2), voxel_coordinates(3)); 
    
    clear AMP_NII
end

%Load in the amplitude data for condition 2
cd(PathName5);
for i = 1:n_group1
    AMP_NII = load_nii(FileName5{1,i});
    AMP_data_cond2(i,:,:,:) = AMP_NII.img;
    
    seedamp_covar(i,2) = AMP_data_cond2(i,voxel_coordinates(1), voxel_coordinates(2), voxel_coordinates(3)); 
    
    clear AMP_NII
end

%Load in the amplitude data for condition 3
cd(PathName6);
for i = 1:n_group1
    AMP_NII = load_nii(FileName6{1,i});
    AMP_data_cond3(i,:,:,:) = AMP_NII.img;
    
    seedamp_covar(i,3) = AMP_data_cond3(i,voxel_coordinates(1), voxel_coordinates(2), voxel_coordinates(3)); 
    
    clear AMP_NII
end



clear i
F_Map_condition = zeros(ref_size1);

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
            
            
            
            
            if sum(COH_data_cond1(:,i,ii,iii)) == 0 || sum(COH_data_cond2(:,i,ii,iii)) == 0 || sum(COH_data_cond3(:,i,ii,iii)) == 0 

                
            elseif all([i,ii,iii] == voxel_coordinates)

                
            else
                %Creating the variables table
                rm_table = table(COH_data_cond1(:,i,ii,iii),...
                                 COH_data_cond2(:,i,ii,iii),...
                                 COH_data_cond3(:,i,ii,iii),...
                                 AMP_data_cond1(:,i,ii,iii),...
                                 AMP_data_cond2(:,i,ii,iii),...
                                 AMP_data_cond3(:,i,ii,iii),...                                 seedamp_covar,...
                                 seedamp_covar(:,1),...
                                 seedamp_covar(:,2),...
                                 seedamp_covar(:,3),...
                                 'VariableNames',...
                                 {'Cond1','Cond2','Cond3','SourceAmp1','SourceAmp2','SourceAmp3','SeedAmp1','SeedAmp2','SeedAmp3'});
            
                rm_design = table([1 2 3]','VariableNames',{'Condition'});
                
                rm_model = ranova(fitrm(rm_table,'Cond1-Cond3~SourceAmp1+SourceAmp2+SourceAmp3+SeedAmp1+SeedAmp2+SeedAmp3','WithinDesign',rm_design));
                
                %The position of the F-value of interest depends on the order that factors were entered into the model
                F_Map_condition(i,ii,iii) = rm_model.F(1);
                
                if stat_counter == 1 %only pull df once
                    
                    df1 = rm_model.DF(1);
                    df2 = rm_model.DF(8);
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


condition_file_name = strcat('ConditionContrast_DF_',df1,'_',df2);



%Save the condition contrast F-map
niftiwrite(F_Map_condition,condition_file_name,nii_info);


total_time = toc;


fprintf('\nThese stats took %4.2f seconds in total\n\n',total_time)

end

