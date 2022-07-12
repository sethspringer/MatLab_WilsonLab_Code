function seth_sourcespace_indAverage_multipleWindows

%PURPOSE:           Create individual averages across several windows (similar to NII Wizard but can do up to six averages)
%
%REQUIRED INPUTS:   NIIs that you want to average.
%

%		  
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, UNMC
%VERSION HISTORY:   03/21/2022  v1: First working version
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Ask the user how many windows they want to average across
inputs = {'How many windows to average?'};
defaults = {'2'};
answer = inputdlg(inputs, 'How many windows/subjects to average?', 2, defaults,'on');
[numb_ave] = deal(answer{:});
numb_ave = str2num(numb_ave);


%Check the inputs to see if they are valid
if isempty(numb_ave)
    error("The input you gave is unacceptable (maybe you entered a letter?)")
elseif isnan(numb_ave)
    error("The input you gave is unacceptable (maybe NaN?)")
elseif ~(numb_ave == 2 | numb_ave == 3 | numb_ave == 4 | numb_ave == 5 | numb_ave == 6)
       error('The input you gave is unacceptable (wrong number probably; must be 2, 3, 4, 5, or 6)')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           LOAD IN DATA                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lood in NIIs and extract information from them
[FileName1,PathName1,~] = uigetfile('*.nii','Select NIIs for Group 1','MultiSelect','on');
cd(PathName1);
% Open first NII and extract parameters%
NII_param1 = load_nii(FileName1{1,1});
if iscell(NII_param1)
    NII_param1 = NII_param1{1,1};
end
ref_size1 = size(NII_param1.img);
nii_info = niftiinfo(FileName1{1,1});


%Lood in NIIs and extract information from them
[FileName2,PathName2,~] = uigetfile('*.nii','Select NIIs for Group 2','MultiSelect','on');
cd(PathName2);
% Open second NII and extract parameters%
NII_param2 = load_nii(FileName2{1,1});
if iscell(NII_param2)
    NII_param2 = NII_param2{1,1};
end
ref_size2 = size(NII_param2.img);


%Lood in NIIs and extract information from them
if numb_ave == 3 | numb_ave == 4 | numb_ave == 5 | numb_ave == 6
    [FileName3,PathName3,~] = uigetfile('*.nii','Select NIIs for Group 3','MultiSelect','on');
    cd(PathName3);
    % Open first NII and extract parameters%
    NII_param3 = load_nii(FileName3{1,1});
    if iscell(NII_param3)
        NII_param3 = NII_param3{1,1};
    end
    ref_size3 = size(NII_param3.img);
    nii_info = niftiinfo(FileName3{1,1});
else
end

%Lood in NIIs and extract information from them
if numb_ave == 4 | numb_ave == 5 | numb_ave == 6
    [FileName4,PathName4,~] = uigetfile('*.nii','Select NIIs for Group 4','MultiSelect','on');
    cd(PathName4);
    % Open first NII and extract parameters%
    NII_param4 = load_nii(FileName4{1,1});
    if iscell(NII_param4)
        NII_param4 = NII_param4{1,1};
    end
    ref_size4 = size(NII_param4.img);
    nii_info = niftiinfo(FileName4{1,1});
else
end


%Lood in NIIs and extract information from them
if numb_ave == 5 | numb_ave == 6
    [FileName5,PathName5,~] = uigetfile('*.nii','Select NIIs for Group 5','MultiSelect','on');
    cd(PathName5);
    % Open first NII and extract parameters%
    NII_param5 = load_nii(FileName5{1,1});
    if iscell(NII_param5)
        NII_param5 = NII_param5{1,1};
    end
    ref_size5 = size(NII_param5.img);
    nii_info = niftiinfo(FileName5{1,1});
else
end

%Lood in NIIs and extract information from them
if numb_ave == 6
    [FileName6,PathName6,~] = uigetfile('*.nii','Select NIIs for Group 6','MultiSelect','on');
    cd(PathName6);
    % Open first NII and extract parameters%
    NII_param6 = load_nii(FileName6{1,1});
    if iscell(NII_param6)
        NII_param6 = NII_param6{1,1};
    end
    ref_size6 = size(NII_param6.img);
    nii_info = niftiinfo(FileName6{1,1});
else
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ENSURE DIMENSIONS FOR ALL NII FILES ARE THE SAME                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check if the dimensions and number of files matches per possible input
if numb_ave == 2
    
    %Check to make sure the NII dimensions and number of subjects is the same%
    if ~isequal(ref_size1,ref_size2)
        error('Volumetric images must be the same dimensions!')
    elseif ~isequal(size(FileName1),size(FileName2))
        error('The same number of files must be loaded!')
    end
    
elseif numb_ave == 3
    
    %Check to make sure the NII dimensions and number of subjects is the same%
    if ~isequal(ref_size1,ref_size2,ref_size3)
        error('Volumetric images must be the same dimensions!')
    elseif ~isequal(size(FileName1),size(FileName2),size(FileName3))
        error('The same number of files must be loaded!')
    end
    
elseif numb_ave == 4
    
    %Check to make sure the NII dimensions and number of subjects is the same%
    if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4)
        error('Volumetric images must be the same dimensions!')
    elseif ~isequal(size(FileName1),size(FileName2),size(FileName3),size(FileName4))
        error('The same number of files must be loaded!')
    end
    
elseif numb_ave == 5
    
    %Check to make sure the NII dimensions and number of subjects is the same%
    if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4,ref_size5)
        error('Volumetric images must be the same dimensions!')
    elseif ~isequal(size(FileName1),size(FileName2),size(FileName3),size(FileName4),size(FileName5))
        error('The same number of files must be loaded!')
    end
    
elseif numb_ave == 6
    
    %Check to make sure the NII dimensions and number of subjects is the same%
    if ~isequal(ref_size1,ref_size2,ref_size3,ref_size4,ref_size5,ref_size6)
        error('Volumetric images must be the same dimensions!')
    elseif ~isequal(size(FileName1),size(FileName2),size(FileName3),size(FileName4),size(FileName5),size(FileName6))
        error('The same number of files must be loaded!')
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                LOAD THE NII DATA INTO 4D MATRICES                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Load the first group of NIIs into 4D matrix (with the first dimension being subject%
cd(PathName1);
for i = 1:size(FileName1,2)
    First_Group_NII = load_nii(FileName1{1,i});
    First_Group_NII_data(i,:,:,:) = First_Group_NII.img;
    clear First_Group_NII
end

clear i %Just to be safe that the next for loop starts fresh

%Load the second group of NIIs into 4D matrix (with the first dimension being subject%
cd(PathName2);
for i = 1:size(FileName2,2)
    Second_Group_NII = load_nii(FileName2{1,i});
    Second_Group_NII_data(i,:,:,:) = Second_Group_NII.img;
    clear Second_Group_NII
end

clear i

if numb_ave == 3 | numb_ave == 4 | numb_ave == 5 | numb_ave == 6
    %Load the third group of NIIs into 4D matrix (with the first dimension being subject%
    cd(PathName3);
    for i = 1:size(FileName3,2)
        Third_Group_NII = load_nii(FileName3{1,i});
        Third_Group_NII_data(i,:,:,:) = Third_Group_NII.img;
        clear Third_Group_NII
    end
else
    
end

clear i

if numb_ave == 4 | numb_ave == 5 | numb_ave == 6
    %Load the fourth group of NIIs into 4D matrix (with the first dimension being subject%
    cd(PathName4);
    for i = 1:size(FileName4,2)
        Fourth_Group_NII = load_nii(FileName4{1,i});
        Fourth_Group_NII_data(i,:,:,:) = Fourth_Group_NII.img;
        clear Fourth_Group_NII
    end
else
    
end

clear i

if numb_ave == 5 | numb_ave == 6
    %Load the fifth group of NIIs into 4D matrix (with the first dimension being subject%
    cd(PathName5);
    for i = 1:size(FileName5,2)
        Fifth_Group_NII = load_nii(FileName5{1,i});
        Fifth_Group_NII_data(i,:,:,:) = Fifth_Group_NII.img;
        clear Fifth_Group_NII
    end
else
    
end

clear i

if numb_ave == 6
    %Load the sixth group of NIIs into 4D matrix (with the first dimension being subject%
    cd(PathName6);
    for i = 1:size(FileName6,2)
        Sixth_Group_NII = load_nii(FileName6{1,i});
        Sixth_Group_NII_data(i,:,:,:) = Sixth_Group_NII.img;
        clear Sixth_Group_NII
    end
else
    
end







save_path = uigetdir(PathName1,'Select the directory to save output average NIIs');
cd(save_path);



NII_Length = (ref_size1(1)*ref_size1(2)*ref_size1(3));
average    = (nan(NII_Length,1))';




%Vectorize the 3D NIIs and average across group
for ii = 1:size(FileName1,2) %Each subject
    
    %Group each subjects' data and remove the extra dimension
    First_Group_NII_data_single  = squeeze(First_Group_NII_data(ii,:,:,:));
    Second_Group_NII_data_single = squeeze(Second_Group_NII_data(ii,:,:,:));
    
    
    
    %Data for groups 3 through 6 are conditionally dealt with
    if numb_ave == 3 | numb_ave == 4 | numb_ave == 5 | numb_ave == 6
        Third_Group_NII_data_single = squeeze(Third_Group_NII_data(ii,:,:,:));
    else
    end
    
    if numb_ave == 4 | numb_ave == 5 | numb_ave == 6
        Fourth_Group_NII_data_single = squeeze(Fourth_Group_NII_data(ii,:,:,:));
    else
    end
    
    if numb_ave == 5 | numb_ave == 6
        Fifth_Group_NII_data_single = squeeze(Fifth_Group_NII_data(ii,:,:,:));
    else
    end
    
    if numb_ave == 6
        Sixth_Group_NII_data_single = squeeze(Sixth_Group_NII_data(ii,:,:,:));
    else
    end
    
    %-------------------------------------------------%
    
    %Perform the average%
    if numb_ave == 2
        for iii = 1:NII_Length %Each position within the 3D
            average(iii) = (First_Group_NII_data_single(iii)+Second_Group_NII_data_single(iii))/2; %Creating an average vector that will need to be reshaped
        end
        
    elseif numb_ave == 3
        for iii = 1:NII_Length %Each position within the 3D
            average(iii) = (First_Group_NII_data_single(iii)+Second_Group_NII_data_single(iii)+Third_Group_NII_data_single(iii))/3; %Creating an average vector that will need to be reshaped
        end
        
    elseif numb_ave == 4
        for iii = 1:NII_Length %Each position within the 3D
            average(iii) = (First_Group_NII_data_single(iii)+Second_Group_NII_data_single(iii)+Third_Group_NII_data_single(iii)+Fourth_Group_NII_data_single(iii))/4; %Creating an average vector that will need to be reshaped
        end
        
    elseif numb_ave == 5
        for iii = 1:NII_Length %Each position within the 3D
            average(iii) = (First_Group_NII_data_single(iii)+Second_Group_NII_data_single(iii)+Third_Group_NII_data_single(iii)+Fourth_Group_NII_data_single(iii)+Fifth_Group_NII_data_single(iii))/5; %Creating an average vector that will need to be reshaped
        end
        
    elseif numb_ave == 6
        for iii = 1:NII_Length %Each position within the 3D
            average(iii) = (First_Group_NII_data_single(iii)+Second_Group_NII_data_single(iii)+Third_Group_NII_data_single(iii)+Fourth_Group_NII_data_single(iii)+Fifth_Group_NII_data_single(iii)+Sixth_Group_NII_data_single(iii))/6; %Creating an average vector that will need to be reshaped
        end
    end
    
    
    
    
    %Reshape the NII into 3Ds of appropriate dimensions
    average_NII = reshape(average,[ref_size1(1),ref_size1(2),ref_size1(3)]);
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  SAVE OUT THE AVERAGE NII FILES                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save_name = strcat(FileName1{ii}(1,1:end-4),'_averaged.nii');
    niftiwrite(average_NII,save_name,nii_info);
    
end %Move on to the next subject



end %end of function


