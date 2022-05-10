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


[FileName1,PathName1,~] = uigetfile('*.nii','Select Coherence NIIs for Group 1','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});

[FileName2,PathName2,~] = uigetfile('*.nii','Select Coherence NIIs for Group 2','MultiSelect','on');
cd(PathName2);
% Open first NII and extract parameters%
NII_param2 = load_nii(FileName2{1,1});
if iscell(NII_param2)
    NII_param2 = NII_param2{1,1};
end
ref_size2 = size(NII_param2.img);












if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4,ref_size5,ref_size6,ref_size7,ref_size8)
    error('Volumetric images must be the same dimensions!')
elseif ~isequal(size(FileName1),size(FileName3)) || ~isequal(size(FileName2),size(FileName4))||~isequal(size(FileName5),size(FileName6))||~isequal(size(FileName7),size(FileName8))
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

[save_name,save_path,save_index] = uiputfile('*.nii','Please save your NII');

tic
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
 
% save_file = NII_param1;
% save_file.img = tMap;
% save_file.hdr.dime.glmax = max(tMap(:));
% save_file.hdr.dime.glmin = min(tMap(:));
cd(save_path);
niftiwrite(tMap,strcat(save_name,save_index),nii_info);
% warning('on','stats:LinearModel:RankDefDesignMat');
total_time = toc;


