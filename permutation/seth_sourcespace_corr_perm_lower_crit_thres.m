function seth_sourcespace_corr_perm_lower_crit_thres




%load the table with all subjects ParID and age...
ages_table = readtable('D:\MIND_VWM\ParIDs_And_Ages.xlsx');



%Run initial correlation to build p-value NII
[FileName1,PathName1,~] = uigetfile('*.nii','Select NIIs to plot pseudo-t values from','MultiSelect','on');
cd(PathName1);

inputs = {'Number of permutations'};
defaults = {'5000'};
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[n_permutations] = deal(answer{:});
n_permutations =str2num(n_permutations);


n_nii = length(FileName1);


%Compute max number of valid permutations (within a reasonable time limit)
%- without resampling!!!
%For a corr or regression, the max number of permutation is factorial(sample size)
%For a paired t-test the max number is factorial(sample size per condition) I believe...

%valid_perms = zeros(n_nii,n_permutations);

for i = 1:n_permutations
    perm_counter = 1;
    while 1
        if perm_counter > n_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',n_permutations,max_perms);
            n_permutations = size(valid_perms,1);
            break
        end
        new_perm = randperm(n_nii);
        perm_counter = perm_counter+1;
        if exist('valid_perms','var') && ~isempty(valid_perms)
            if ~ismember(new_perm,valid_perms,'rows') & ~isequal(new_perm,[1:n_nii])
                valid_perms(i,:) = new_perm;
                break
            end
        elseif ~isequal(new_perm,[1:n_nii])
            valid_perms(i,:) = new_perm;
        end
    end
end



% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});


%Load in group 1 data (controls)
cd(PathName1);
for i = 1:n_nii
    nii_data_temp = load_nii(FileName1{1,i});
    nii_data(i,:,:,:) = nii_data_temp.img;
    clear average_nii_temp
end


ages = ages_table.Age;
r_map = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);
p_map = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);
sig_masking_map = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);

for x = 1:ref_size1(1)
    for y = 1:ref_size1(2)
        for z = 1:ref_size1(3)
            
            if sum(nii_data(:,x,y,z)) ~= 0
                
                nii_data_vector = nii_data(:,x,y,z);
                
                [r,p] = corr(nii_data_vector,ages);
                
                r_map(x,y,z) = r;
                p_map(x,y,z) = p;
                
                
            end
            
        end
    end
end

%Create a map where you have significant voxels marked
sig_masking_map(p_map <= 0.025 & p_map > 0 & r_map ~= 0) = 1;

if nnz(sig_masking_map) == 0
    error('There are no significant peaks')
end

sig_masking_map_cluster_building = sig_masking_map; %copy this so you can change the values when cluster building


cluster_label = 1000; %This will be used later. All clusters will be labeled in the thousands...


clear i

%Scan through (in the x plane) the sig brain for 1s, when you find one stop and build a cluster
for x_plane = 1:ref_size1(1)
    
    current_x_plane = squeeze(sig_masking_map_cluster_building(x_plane,:,:));
    
    if sum(current_x_plane,'all') > 0 %Meaning you have found a significant cluster
        
        [y_points,z_points] = find(current_x_plane==1); %Find the location of all 1s in this x_plane
        
        
        for i = 1:length(y_points) %now you need to build a cluster for each significant voxel you found in this slice
            
            if sig_masking_map_cluster_building(x_plane,y_points(i),z_points(i)) == 1 %This is making sure that the voxel you are on is not part of another cluster...
                
                %%%%%%% Build the cluster(s)%%%%%%%%
                %%%%%%% As the clusters grow, each new layer gets assigned a new value%%%%
                
                sig_masking_map_cluster_building(x_plane,y_points(i),z_points(i)) = 2; %Change the starting voxel to 2, just so it is not 1 anymore
                
                cluster_iteration = 0;
                
                end_of_cluster_building = 1; %needs to not equal 0 at first so the while loop starts
                
                while end_of_cluster_building ~= 0 %
                    
                    end_of_cluster_building = 0; %Always reset this to zero so you break the while loop if the cluster doesn't grow
                    cluster_iteration = cluster_iteration+1; %Keep track of the number of cluster layers
                    
                    
                    %Find values for the last round of cluster building (if this is the first round, find the starting voxel)
                    %This is the way to do a 3D find matrix
                    clear size
                    f = find (sig_masking_map_cluster_building == cluster_iteration+1);
                    [x_for_this_cluster_round y_for_this_cluster_round z_for_this_cluster_round] = ind2sub(size(sig_masking_map_cluster_building),f);
                    
                    %This is a list of all of the voxel coordinates that you have to work from this cluster round
                    coords_for_this_round = [x_for_this_cluster_round,y_for_this_cluster_round,z_for_this_cluster_round];
                    
                    %need to loop over each of these voxels that were found in the last round and find their significant neighbors
                    for ii = 1:size(coords_for_this_round,1)
                        
                            %%%%%%%%%% First check the voxels that share a side
                            %+x
                            if ~(coords_for_this_round(ii,1)+1 > ref_size1(1)) %Only check +x if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x
                            if ~(coords_for_this_round(ii,1)-1 == 0) %Only check -x if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %+y
                            if ~(coords_for_this_round(ii,2)+1 > ref_size1(2)) %Only check +y if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-y
                            if ~(coords_for_this_round(ii,2)-1 == 0) %Only check -y if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %+z
                            if ~(coords_for_this_round(ii,3)+1 > ref_size1(3)) %Only check +y if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-z
                            if ~(coords_for_this_round(ii,3)-1 == 0) %Only check -z if it doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            %%%%%%%%%% Second check the voxels that share an edge
                            %+x/+y
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,2)+1 > ref_size1(2))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+x/-y
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,2)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-x/+y
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,2)+1 > ref_size1(2))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x/-y
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,2)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            
                            %+x/+z
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+x/-z
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-x/+z
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x/-z
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            
                            %+y/+z
                            if ~((coords_for_this_round(ii,2)+1 > ref_size1(2)) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+y/-z
                            if ~((coords_for_this_round(ii,2)+1 > ref_size1(2)) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-y/+z
                            if ~((coords_for_this_round(ii,2)-1 == 0) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-y/-z
                            if ~((coords_for_this_round(ii,2)-1 == 0) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    end %end of the cluster-building for loop
                    
                end %end of cluster-building while loop
                
                
                
                
            end %End of if statement that starts the cluster building by making sure the starting voxel has not been part of another cluster
            
            %Once you have built a cluster, you should label as unique before moving on to builind the next one
            cluster_label = cluster_label + 1;
            
            sig_masking_map_cluster_building(sig_masking_map_cluster_building>1.5 & sig_masking_map_cluster_building < 999) = cluster_label;
            
            
            
        end
        
    end
    
end %end of the looping through the x-planes




%     niftiwrite(sig_masking_map_cluster_building,'test_clustering_n4400_n4000.nii',nii_info)
%
%     unique(sig_masking_map_cluster_building)





%Now that you have the clusters separated, you need to sum the test stats for each cluster
cluster_labels = unique(sig_masking_map_cluster_building(sig_masking_map_cluster_building~=0)); %list the cluster labels
n_clusters = length(cluster_labels);

if n_clusters == 0
    error('There are no significant clusters, try a less stringent p-value threshold')
end

clear i
for i = 1:n_clusters
    
    cluster_sum(i) = sum(r_map(sig_masking_map_cluster_building == cluster_labels(i))); %Extract the test stats for each cluster
    
end



























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                                 PERMUTATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cluster_sum_perm_max_list = zeros(n_permutations,1);


%Initialize the loading bar
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0);
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);

tic

%Now I can move on to permuting and sampling clusters...
for permutation_index = 1:n_permutations
    
    
    r_map_perm = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);
    p_map_perm = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);
    sig_masking_map_perm = zeros([ref_size1(1),ref_size1(2),ref_size1(3)]);
    
    
    
    %I am shuffling the ages rather than the NIIs, though the outcome is the same
    ages_perm = ages(valid_perms(permutation_index,:));
    
    
    
    for x = 1:ref_size1(1)
        for y = 1:ref_size1(2)
            for z = 1:ref_size1(3)
                
                if sum(nii_data(:,x,y,z)) ~= 0
                    
                    nii_data_vector = nii_data(:,x,y,z);
                    
                    [r,p] = corr(nii_data_vector,ages_perm);
                    
                    r_map_perm(x,y,z) = r;
                    p_map_perm(x,y,z) = p;
                    
                    
                end
                
            end
        end
    end
    
    %Create a map where you have significant voxels marked
    sig_masking_map_perm(p_map_perm <= 0.025 & p_map_perm > 0 & r_map_perm ~= 0) = 1;
    
    if nnz(sig_masking_map_perm) == 0
        continue %if there are no significant voxels for this permutation, go on to the next one.
    end
    
    sig_masking_map_perm_cluster_building = sig_masking_map_perm; %copy this so you can change the values when cluster building
    
    
    cluster_label_perm = 1000; %This will be used later. All clusters will be labeled in the thousands...
    
    
    clear i
    
    %Scan through (in the x plane) the sig brain for 1s, when you find one stop and build a cluster
    
    for x_plane = 1:ref_size1(1)
        
        current_x_plane = squeeze(sig_masking_map_perm_cluster_building(x_plane,:,:));
        
        if sum(current_x_plane,'all') > 0 %Meaning you have found a significant cluster
            
            [y_points,z_points] = find(current_x_plane==1); %Find the location of all 1s in this x_plane
            
            
            for i = 1:length(y_points) %now you need to build a cluster for each significant voxel you found in this slice
                
                if sig_masking_map_perm_cluster_building(x_plane,y_points(i),z_points(i)) == 1 %This is making sure that the voxel you are on is not part of another cluster...
                    
                    %%%%%%% Build the cluster(s)%%%%%%%%
                    %%%%%%% As the clusters grow, each new layer gets assigned a new value%%%%
                    
                    sig_masking_map_perm_cluster_building(x_plane,y_points(i),z_points(i)) = 2; %Change the starting voxel to 2, just so it is not 1 anymore
                    
                    cluster_iteration = 0;
                    
                    end_of_cluster_building = 1; %needs to not equal 0 at first so the while loop starts
                    
                    while end_of_cluster_building ~= 0 %
                        
                        end_of_cluster_building = 0; %Always reset this to zero so you break the while loop if the cluster doesn't grow
                        cluster_iteration = cluster_iteration+1; %Keep track of the number of cluster layers
                        
                        
                        %Find values for the last round of cluster building (if this is the first round, find the starting voxel)
                        %This is the way to do a 3D find matrix
                        clear size
                        f = find (sig_masking_map_perm_cluster_building == cluster_iteration+1);
                        [x_for_this_cluster_round y_for_this_cluster_round z_for_this_cluster_round] = ind2sub(size(sig_masking_map_perm_cluster_building),f);
                        
                        %This is a list of all of the voxel coordinates that you have to work from this cluster round
                        coords_for_this_round = [x_for_this_cluster_round,y_for_this_cluster_round,z_for_this_cluster_round];
                        
                        %need to loop over each of these voxels that were found in the last round and find their significant neighbors
                        for ii = 1:size(coords_for_this_round,1)
                            
                            %%%%%%%%%% First check the voxels that share a side
                            %+x
                            if ~(coords_for_this_round(ii,1)+1 > ref_size1(1)) %Only check +x if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x
                            if ~(coords_for_this_round(ii,1)-1 == 0) %Only check -x if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %+y
                            if ~(coords_for_this_round(ii,2)+1 > ref_size1(2)) %Only check +y if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-y
                            if ~(coords_for_this_round(ii,2)-1 == 0) %Only check -y if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %+z
                            if ~(coords_for_this_round(ii,3)+1 > ref_size1(3)) %Only check +y if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-z
                            if ~(coords_for_this_round(ii,3)-1 == 0) %Only check -z if it doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            %%%%%%%%%% Second check the voxels that share an edge
                            %+x/+y
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,2)+1 > ref_size1(2))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+x/-y
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,2)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-x/+y
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,2)+1 > ref_size1(2))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x/-y
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,2)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            
                            %+x/+z
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+x/-z
                            if ~((coords_for_this_round(ii,1)+1 > ref_size1(1)) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)+1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-x/+z
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-x/-z
                            if ~((coords_for_this_round(ii,1)-1 == 0) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1)-1,coords_for_this_round(ii,2),coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            
                            %+y/+z
                            if ~((coords_for_this_round(ii,2)+1 > ref_size1(2)) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %+y/-z
                            if ~((coords_for_this_round(ii,2)+1 > ref_size1(2)) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)+1,coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            %-y/+z
                            if ~((coords_for_this_round(ii,2)-1 == 0) || (coords_for_this_round(ii,3)+1 > ref_size1(3))) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)+1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)+1) = cluster_iteration+2; %each voxel added to the cluster should increase in value
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1; %change this from zero, showing the cluster is still growing
                                end
                            end
                            %-y/-z
                            if ~((coords_for_this_round(ii,2)-1 == 0) || (coords_for_this_round(ii,3)-1 == 0)) %Only check +x and +y if they doesn't index out of bounds
                                if sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)-1) == 1
                                    sig_masking_map_perm_cluster_building(coords_for_this_round(ii,1),coords_for_this_round(ii,2)-1,coords_for_this_round(ii,3)-1) = cluster_iteration+2;
                                    
                                    end_of_cluster_building = end_of_cluster_building + 1;
                                end
                            end
                            
                            
                            
                            
                            
                            
                            
                            
                            
                        end %end of the cluster-building for loop
                        
                    end %end of cluster-building while loop
                    
                    
                    
                    
                end %End of if statement that starts the cluster building by making sure the starting voxel has not been part of another cluster
                
                %Once you have built a cluster, you should label as unique before moving on to builind the next one
                cluster_label_perm = cluster_label_perm + 1;
                
                sig_masking_map_perm_cluster_building(sig_masking_map_perm_cluster_building>1.5 & sig_masking_map_perm_cluster_building < 999) = cluster_label_perm;
                
                
                
            end
            
        end
        
    end %end of the looping through the x-planes
    
    
    
    
    %     niftiwrite(sig_masking_map_perm_cluster_building,'test_clustering_new.nii',nii_info)
    %
    %     unique(sig_masking_map_perm_cluster_building)
    
    clear n_clusters
    
    %Now that you have the clusters separated, you need to sum the test stats for each cluster
    cluster_labels_perm = unique(sig_masking_map_perm_cluster_building(sig_masking_map_perm_cluster_building~=0)); %list the cluster labels
    n_clusters = length(cluster_labels_perm);
    
    
    cluster_sum_perm = zeros(n_clusters,1);
    
    if n_clusters > 0
        clear i
        for i = 1:n_clusters
            
            cluster_sum_perm(i) = sum(r_map_perm(sig_masking_map_perm_cluster_building == cluster_labels_perm(i))); %Extract the test stats for each cluster
            
        end
    end
    
    %Now you have all of the clusters for this permutation, find the absolute max, and move on to the next permutation...
    if n_clusters == 0
        cluster_sum_perm_max_list(permutation_index) = 0; %This is redundant
    else
        [~,cluster_sum_perm_max_index] = max(abs(cluster_sum_perm)); %find the position of abs(max())
        cluster_sum_perm_max_list(permutation_index) = cluster_sum_perm(cluster_sum_perm_max_index); %grab the abs(max())
    end
    
    
if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Permutation testing canceled by user')
end
waitbar(permutation_index/n_permutations,wait_bar,sprintf('Permutation #%.0f',permutation_index))

end %end of permutation for loop

delete(wait_bar);

total_time = toc



%Now you have the null distribution from the permutaiton testing, compare the clusters from the original data...
%%%%% This is how Alex did it in the timeseries Monte Carlo toolbox, this
%%%%% needs to be futher investigated.....

cluster_sum_perm_max_list_positive = (cluster_sum_perm_max_list(cluster_sum_perm_max_list >= 0));
cluster_sum_perm_max_list_negative = (cluster_sum_perm_max_list(cluster_sum_perm_max_list <= 0));

%calculate the p-value of each cluster from the original timeseries
cluster_pvals = zeros(1,size(cluster_sum,2));
for i = 1:size(cluster_sum,2) 
    if cluster_sum(1,i)>0
        cluster_pvals(1,i) = 1-((sum(sum(cluster_sum_perm_max_list_positive<cluster_sum(1,i)))/sum(sum(~isnan(cluster_sum_perm_max_list_positive)))));
    else
        cluster_pvals(1,i) = 1-((sum(sum(cluster_sum_perm_max_list_negative>cluster_sum(1,i)))/sum(sum(~isnan(cluster_sum_perm_max_list_negative)))));
    end
end



%Now you have the null distribution from the permutaiton testing, compare the clusters from the original data...
%%%%% This is using both side of the distribution...

%calculate the p-value of each cluster from the original timeseries
cluster_pvals_test_new = zeros(1,size(cluster_sum,2));
for i = 1:size(cluster_sum,2) 
    if cluster_sum(1,i)>0
        cluster_pvals_test_new(1,i) = sum(cluster_sum_perm_max_list>cluster_sum(1,i))/sum(~isnan(cluster_sum_perm_max_list));
    else
        cluster_pvals_test_new(1,i) = sum(cluster_sum_perm_max_list<cluster_sum(1,i))/sum(~isnan(cluster_sum_perm_max_list));
    end
end




%Another way, from https://www.biomedware.com/files/documentation/OldCSHelp/MCR/Calculating_Monte_Carlo_p-values.htm

cluster_pvals_test_new_new = zeros(1,size(cluster_sum,2));
for i = 1:size(cluster_sum,2) 
    if cluster_sum(1,i)>0
        cluster_pvals_test_new_new(1,i) = (sum(cluster_sum_perm_max_list>cluster_sum(1,i))+1)/(sum(~isnan(cluster_sum_perm_max_list))+1);
    else
        cluster_pvals_test_new_new(1,i) = (sum(cluster_sum_perm_max_list<cluster_sum(1,i))+1)/(sum(~isnan(cluster_sum_perm_max_list))+1);
    end
end






%Now create outputs and sort them
cluster_sum = cluster_sum';
cluster_pvals_test_new_new = cluster_pvals_test_new_new';

%Need to decide with p-values to use... xx
sorting_list = array2table([cluster_labels,cluster_sum,cluster_pvals_test_new_new]);

sorting_list_sorted = sortrows(sorting_list,'Var3');

sorted_p_values       = sorting_list_sorted.Var3;
sorted_cluster_size   = sorting_list_sorted.Var2;
sorted_cluster_labels = sorting_list_sorted.Var1;

new_cluster_labels    = (1:length(sorted_cluster_labels))';

output_table = array2table([new_cluster_labels,sorted_cluster_size,sorted_p_values]);


%Relabel the output mask values to matach the output table
output_nii = sig_masking_map_cluster_building;

clear i
for i = 1:length(sorted_cluster_labels)
    
    output_nii(output_nii == sorted_cluster_labels(i)) = new_cluster_labels(i);
    
end




cd(PathName1);
niftiwrite(output_nii,'cluster_permutation_mask.nii',nii_info)









end %end of function