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


%Preallocate
COH_data_cond1 = zeros([(sum(num_per_group)),ref_size(1,:)]);
COH_data_cond2 = zeros([(sum(num_per_group)),ref_size(1,:)]);

AMP_data_cond1 = zeros([(sum(num_per_group)),ref_size(1,:)]);
AMP_data_cond2 = zeros([(sum(num_per_group)),ref_size(1,:)]);


pred_groups = [zeros(num_per_group(1),1);ones(num_per_group(2),1);ones(num_per_group(3),1)*-1;ones(num_per_group(4),1)*-2];
seedamp_covar_cond1 = zeros(sum(num_per_group),1);
seedamp_covar_cond2 = zeros(sum(num_per_group),1);


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





tMap = zeros(ref_size1);
% warning('off','stats:LinearModel:RankDefDesignMat');
for i = 1:ref_size1(1,1)
    fprintf('Running X:%d\n',i)
    for ii = 1:ref_size1(1,2)
        for iii = 1:ref_size1(1,3)
            dep = COH_data(:,i,ii,iii);
            voxelamp_covar = AMP_data(:,i,ii,iii);
            if sum(dep(:)) == 0 || sum(voxelamp_covar(:)) == 0
                tMap(i,ii,iii) = 0;
            elseif all([i,ii,iii] == voxel_coordinates)
                tMap(i,ii,iii) = 0;
            else
                pred = [pred_groups,voxelamp_covar,seedamp_covar];
                model = fitglm(pred,dep);
                tMap(i,ii,iii) = model.Coefficients.tStat(2,1);
            end
        end
    end
end

% save_file = NII_param;
% save_file.img = tMap;
% save_file.hdr.dime.glmax = max(tMap(:));
% save_file.hdr.dime.glmin = min(tMap(:));
cd(save_path);
niftiwrite(tMap,strcat(save_name,save_index),nii_info);
% warning('on','stats:LinearModel:RankDefDesignMat');
total_time = toc;


