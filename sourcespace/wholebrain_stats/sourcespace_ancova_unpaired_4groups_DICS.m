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

%Pick up here...
COH_data = zeros([(size(FileName1,2)+size(FileName2,2)+size(FileName5,2)+size(FileName7,2)),ref_size1]);
AMP_data = zeros([(size(FileName3,2)+size(FileName4,2)+size(FileName6,2)+size(FileName8,2)),ref_size1]);
% COH_data = [zeros([size(FileName1,2),ref_size1]);zeros([size(FileName2,2),ref_size2])];
% AMP_data = [zeros([size(FileName3,2),ref_size3]);zeros([size(FileName4,2),ref_size4])];
pred_groups = [zeros(size(FileName1,2),1);ones(size(FileName2,2),1);ones(size(FileName5,2),1)*-1;ones(size(FileName7,2),1)*-2];
seedamp_covar = zeros(size(FileName1,2)+size(FileName2,2)+size(FileName5,2)+size(FileName7,2),1);



cd(PathName1);
for i = 1:size(FileName1,2)
    COH_NII = load_nii(FileName1{1,i});
    COH_data(i,:,:,:) = COH_NII.img;
    clear COH_NII
end

counter = 1;
cd(PathName2);
for i = (size(FileName1,2)+1):(size(FileName1,2)+size(FileName2,2))
    COH_NII = load_nii(FileName2{1,counter});
    COH_data(i,:,:,:) = COH_NII.img;
    clear COH_NII
    counter = counter+1;
end


counter = 1;
cd(PathName5);
for i = (size(FileName1,2)+size(FileName2,2)+1):(size(FileName1,2)+size(FileName2,2)+size(FileName5,2))
    COH_NII = load_nii(FileName5{1,counter});
    COH_data(i,:,:,:) = COH_NII.img;
    clear COH_NII
    counter = counter+1;
end
counter = 1;
cd(PathName7);
for i = (size(FileName1,2)+size(FileName2,2)+size(FileName5,2)+1):(size(FileName1,2)+size(FileName2,2)+size(FileName5,2)+size(FileName7,2))
    COH_NII = load_nii(FileName7{1,counter});
    COH_data(i,:,:,:) = COH_NII.img;
    clear COH_NII
    counter = counter+1;
end

cd(PathName3);
for i = 1:size(FileName3,2)
    AMP_NII = load_nii(FileName3{1,i});
    AMP_data(i,:,:,:) = AMP_NII.img;
    seedamp_covar(i,1) = AMP_data(i,voxel_coordinates(1,1), voxel_coordinates(1,2), voxel_coordinates(1,3));
    clear AMP_NII
end

counter = 1;
cd(PathName4);
for i = (size(FileName3,2)+1):(size(FileName3,2)+size(FileName4,2))
    AMP_NII = load_nii(FileName4{1,counter});
    AMP_data(i,:,:,:) = AMP_NII.img;
    seedamp_covar(i,1) = AMP_data(i,voxel_coordinates(1,1), voxel_coordinates(1,2), voxel_coordinates(1,3));
    clear AMP_NII
    counter = counter+1;
end


counter = 1;
cd(PathName6);
for i = (size(FileName3,2)+size(FileName4,2)+1):(size(FileName3,2)+size(FileName4,2)+size(FileName6,2))
    AMP_NII = load_nii(FileName6{1,counter});
    AMP_data(i,:,:,:) = AMP_NII.img;
    seedamp_covar(i,1) = AMP_data(i,voxel_coordinates(1,1), voxel_coordinates(1,2), voxel_coordinates(1,3));
    clear AMP_NII
    counter = counter+1;
end

counter = 1;
cd(PathName8);
for i = (size(FileName3,2)+size(FileName4,2)+1)+size(FileName6,2):(size(FileName3,2)+size(FileName4,2)+size(FileName6,2)+size(FileName8,2))
    AMP_NII = load_nii(FileName8{1,counter});
    AMP_data(i,:,:,:) = AMP_NII.img;
    seedamp_covar(i,1) = AMP_data(i,voxel_coordinates(1,1), voxel_coordinates(1,2), voxel_coordinates(1,3));
    clear AMP_NII
    counter = counter+1;
end



% seedamp_covar = [seedamp1,seedamp2]';

% if ~exist('seedamp_covar','var')
%     [SeedAmp1File,SeedAmp1Path,~] = uigetfile('*.txt','Select File for Seed Amp for Group 1','MultiSelect','on');
%     [SeedAmp2File,SeedAmp2Path,~] = uigetfile('*.txt','Select File for Seed Amp for Group 2','MultiSelect','on');
%     cd(SeedAmp1Path);
%     fid = fopen(SeedAmp1File,'rt');
%     seedamp1 = cell2mat(textscan(fid, '%*s %f %*[^\n]','HeaderLines',1));
%     fclose(fid);
%     cd(SeedAmp2Path);
%     fid = fopen(SeedAmp2File,'rt');
%     seedamp2 = cell2mat(textscan(fid, '%*s %f %*[^\n]','HeaderLines',1));
%     fclose(fid);
%     seedamp_covar = [seedamp1;seedamp2];
% end
% 
% if ~exist('SeedAmp1File','var') || ~exist('SeedAmp2File','var')
%     input_titles = {'Seed AMP G1','Seed AMP G2'};
%     default = {'',''};
%     answer = inputdlg(input_titles, 'Please Input Seed AMP Data', 2, default,'on');
%     seedamp_covar = [str2num(answer{1,1});str2num(answer{2,1})];
% end

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


