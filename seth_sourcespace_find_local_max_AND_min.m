function seth_sourcespace_find_local_max_AND_min

%PURPOSE:           Find pseudo-t values for local maxima and minima NIIs.
%
%REQUIRED INPUTS:   Average NIIs.
%
%		    
%
%NOTES:            
%		  
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   03/31/2022  v1: First working version of program




[FileName1,PathName1,~] = uigetfile('*.nii','Select Average NIIs','MultiSelect','on');
cd(PathName1);



if ~iscell(FileName1)
    files = {FileName1};
end

% Open first NII and extract parameters%
NII = load_nii(FileName1);
NII_data = NII.img;
ref_size1 = size(NII.img);
nii_info = niftiinfo(FileName1);






%%%%%%%%%%%%%%%%% Indexing coordinates to check all surrounding positions %%%%%%%%%%%%%%%%%
%
% %6 sides
% %12 edges
% %8 corners
%
% %SIDES
%
% %x plane
% NII_data(i+1,ii,iii)
% NII_data(i-1,ii,iii)
%
% %y plane
% NII_data(i,ii+1,iii)
% NII_data(i,ii-1,iii)
%
% %z plane
% NII_data(i,ii,iii+1)
% NII_data(i,ii,iii-1)
%
%
% %EDGES
%
% %x/y plane
% NII_data(i+1,ii+1,iii)
% NII_data(i+1,ii-1,iii)
% NII_data(i-1,ii+1,iii)
% NII_data(i-1,ii-1,iii)
%
% %x/z plane
% NII_data(i+1,ii,iii+1)
% NII_data(i-1,ii,iii+1)
% NII_data(i-1,ii,iii-1)
% NII_data(i+1,ii,iii-1)
%
% %y/z plane
% NII_data(i,ii+1,iii+1)
% NII_data(i,ii-1,iii+1)
% NII_data(i,ii-1,iii-1)
% NII_data(i,ii+1,iii-1)
%
%
%
% %Corners - In testing, I found that neuroelf does not consider voxels
% that only touch corners to be neighbors...
%
% %+z
% NII_data(i+1,ii+1,iii+1)
% NII_data(i-1,ii+1,iii+1)
% NII_data(i-1,ii-1,iii+1)
% NII_data(i+1,ii-1,iii+1)
%
% %-z
% NII_data(i+1,ii+1,iii-1)
% NII_data(i-1,ii+1,iii-1)
% NII_data(i-1,ii-1,iii-1)
% NII_data(i+1,ii-1,iii-1)
%

%I have noticed an issue where at the edge of the brain are messing with
%the conditionals in the folloding loop...
NII_data(isnan(NII_data))=0;


pseudo_t_values_of_max_peaks = 0;
max_peak_index = 1;

%For each x, y, and z you have to start and end one away from the edge so
%the index doesn't exceed array bounds
%FINDING LOCAL MAXIMA
for i = 2:ref_size1(1,1)-1
    for ii = 2:ref_size1(1,2)-1
        for iii = 2:ref_size1(1,3)-1
            
            if NII_data(i,ii,iii) == 0
                %Do nothing
            else
                
                %Check voxels that share a side
                if (NII_data(i,ii,iii) > NII_data(i+1,ii,iii)) && (NII_data(i,ii,iii) > NII_data(i-1,ii,iii)) &&... %Check the x plane
                        (NII_data(i,ii,iii) > NII_data(i,ii+1,iii)) && (NII_data(i,ii,iii) > NII_data(i,ii-1,iii)) &&... %Check the y plane
                        (NII_data(i,ii,iii) > NII_data(i,ii,iii+1)) && (NII_data(i,ii,iii) > NII_data(i,ii,iii-1))       %Check the z plane
                    
                    %Only continues to the voxels that share an edge if the side check is cleared
                    if (NII_data(i,ii,iii) > NII_data(i+1,ii+1,iii)) && (NII_data(i,ii,iii) > NII_data(i+1,ii-1,iii)) && (NII_data(i,ii,iii) > NII_data(i-1,ii+1,iii)) && (NII_data(i,ii,iii) > NII_data(i-1,ii-1,iii)) &&...    %x/y plane
                            (NII_data(i,ii,iii) > NII_data(i+1,ii,iii+1)) && (NII_data(i,ii,iii) > NII_data(i-1,ii,iii+1)) && (NII_data(i,ii,iii) > NII_data(i-1,ii,iii-1)) && (NII_data(i,ii,iii) > NII_data(i+1,ii,iii-1)) &&...    %x/z plane
                            (NII_data(i,ii,iii) > NII_data(i,ii+1,iii+1)) && (NII_data(i,ii,iii) > NII_data(i,ii-1,iii+1)) && (NII_data(i,ii,iii) > NII_data(i,ii-1,iii-1)) && (NII_data(i,ii,iii) > NII_data(i,ii+1,iii-1))          %y/z plane
                        
                        
                        pseudo_t_values_of_max_peaks(max_peak_index,1) = NII_data(i,ii,iii);
                        max_peak_index = max_peak_index+1;
                        
                    else
                        %Do nothing
                    end
                else
                    %Do nothing
                end
            end
        end
    end
end


ordered_max_pseudo_ts = unique(pseudo_t_values_of_max_peaks);

max_filename = strcat(FileName1(1:end-4),'_local_maximum_values.txt');

dlmwrite(max_filename,ordered_max_pseudo_ts);










clear i
clear ii
clear iii



pseudo_t_values_of_min_peaks = 0;
min_peak_index = 1;

%For each x, y, and z you have to start and end one away from the edge so
%the index doesn't exceed array bounds
%FINDING LOCAL MINIMA
for i =2:ref_size1(1,1)-1
    for ii = 2:ref_size1(1,2)-1
        for iii = 2:ref_size1(1,3)-1
            
            if NII_data(i,ii,iii) == 0
                %Do nothing
            else
                
                %Check voxels that share a side
                if (NII_data(i,ii,iii) < NII_data(i+1,ii,iii)) && (NII_data(i,ii,iii) < NII_data(i-1,ii,iii)) &&... %Check the x plane
                        (NII_data(i,ii,iii) < NII_data(i,ii+1,iii)) && (NII_data(i,ii,iii) < NII_data(i,ii-1,iii)) &&... %Check the y plane
                        (NII_data(i,ii,iii) < NII_data(i,ii,iii+1)) && (NII_data(i,ii,iii) < NII_data(i,ii,iii-1))       %Check the z plane
                    
                    %Only continues to the voxels that share an edge if the side check is cleared
                    if (NII_data(i,ii,iii) < NII_data(i+1,ii+1,iii)) && (NII_data(i,ii,iii) < NII_data(i+1,ii-1,iii)) && (NII_data(i,ii,iii) < NII_data(i-1,ii+1,iii)) && (NII_data(i,ii,iii) < NII_data(i-1,ii-1,iii)) &&...    %x/y plane
                            (NII_data(i,ii,iii) < NII_data(i+1,ii,iii+1)) && (NII_data(i,ii,iii) < NII_data(i-1,ii,iii+1)) && (NII_data(i,ii,iii) < NII_data(i-1,ii,iii-1)) && (NII_data(i,ii,iii) < NII_data(i+1,ii,iii-1)) &&...    %x/z plane
                            (NII_data(i,ii,iii) < NII_data(i,ii+1,iii+1)) && (NII_data(i,ii,iii) < NII_data(i,ii-1,iii+1)) && (NII_data(i,ii,iii) < NII_data(i,ii-1,iii-1)) && (NII_data(i,ii,iii) < NII_data(i,ii+1,iii-1))          %y/z plane
                        
                        
                        pseudo_t_values_of_min_peaks(min_peak_index,1) = NII_data(i,ii,iii);
                        min_peak_index = min_peak_index+1;
                        
                    else
                        %Do nothing
                    end
                else
                    %Do nothing
                end
            end
        end
    end
end


ordered_min_pseudo_ts = unique(pseudo_t_values_of_min_peaks);

min_filename = strcat(FileName1(1:end-4),'_local_minimum_values.txt');

dlmwrite(min_filename,ordered_min_pseudo_ts);


