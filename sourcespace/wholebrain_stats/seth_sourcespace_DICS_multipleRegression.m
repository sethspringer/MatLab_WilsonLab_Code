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
%VERSION HISTORY:   12/1/2022  v1: First working version of program


%load the table with all subjects ParID and age...
ages_table = readtable('D:\MIND_VWM\ParIDs_And_Ages.xlsx');


[FileName1,PathName1,~] = uigetfile('*.nii','Select DICS NIIs','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});

[FileName2,PathName2,~] = uigetfile('*.nii','Select AMP NIIs','MultiSelect','on');
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


[save_name,save_path] = uiputfile('.nii','Select the directory to save the statistical map and give a name');



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






n_nii = size(FileName1,2);


%Preallocate nii and seed amp arrayes
COH_data_cond1 = zeros([n_nii,ref_size1]);
COH_data_cond2 = zeros([n_nii,ref_size1]);
seedamp_covar = zeros(n_nii,1);



%Load in DICS NIIs
cd(PathName1);
for i = 1:n_nii
    COH_NII = load_nii(FileName1{1,i});
    COH_data_cond1(i,:,:,:) = COH_NII.img;
    clear COH_NII
end


clear i

%Load in AMP NIIs; also grab the AMP values from the seed voxel
cd(PathName2);
for i = 1:n_nii
    COH_NII = load_nii(FileName2{1,i});
    COH_data_cond2(i,:,:,:) = COH_NII.img;
    seedamp_covar(i,1) = COH_data_cond2(i,voxel_coordinates(1,1), voxel_coordinates(1,2), voxel_coordinates(1,3));
    
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
            if sum(COH_data_cond1(:,i,ii,iii)) == 0 || sum(COH_data_cond2(:,i,ii,iii)) == 0 %can't run stats on zeros
                
            elseif all([i,ii,iii] == voxel_coordinates) %dont run stats on seed
                
            else
                
                
                
                
                
                %Creating the variables table
                variable_table = table(...
                    COH_data_cond1(:,i,ii,iii),...
                    COH_data_cond2(:,i,ii,iii),...
                    seedamp_covar,...
                    ages_OE,...
                    'VariableNames',...
                    {'DICS','AMP','seed_amp','age'});
                
                model = fitlm(variable_table, 'DICS~age+AMP+seed_amp');
                
               

                %The position of the F-value of interest depends on the order that factors were entered into the model
                F_Map_condition(i,ii,iii) =  model.Coefficients.tStat(4);
                
                if stat_counter == 1; %only pull df once
                    
                    df = model.DFE;
                    stat_counter = stat_counter + 1; %Use this to print out the df text file
                    
                end
                
                
            end
        end
    end
end

cd(save_path);

df = num2str(df);


close(progress_bar)


save_name = strcat(save_name(1:end-4),'_DF_',df);

%Save the condition contrast F-map
niftiwrite(F_Map_condition,save_name,nii_info);

total_time = toc;


fprintf('\nThese stats took %4.2f seconds in total\n\n',total_time)

end

