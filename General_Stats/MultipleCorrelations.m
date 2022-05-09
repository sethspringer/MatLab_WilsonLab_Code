function MultipleCorrelations(first_x_column,last_x_column,first_y_column,last_y_column)

%PURPOSE:          Correlate one group of variables with another.
%
%REQUIRED INPUTS:  An excel spreadsheet organized like it would be in SPSS.
%                  You MUST specify what columns that each variable begins
%                  and ends in...
%
%NOTES:            Removes missing values pair-wise.
%                  
%                  
%AUTHOR:           Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:  08/31/2020  v1: First working version of program
%



[file,path] = uigetfile('*.xlsx','Select the excel you want to work with...');

cd(path)

%Read in the excel file
[numData,txtData,rawData] = xlsread(file);


%Pull the x variable headers
x_variable_headers = rawData(1,first_x_column:last_x_column);
%Pull the y variable headers
y_variable_headers = rawData(1,first_y_column:last_y_column)';

%Pull the data for the x variables
x_variables         = numData(:,first_x_column:last_x_column);
%Determine how many x variables there are
number_x_variables  = size(x_variables,2);
%Determine the AVE for each x variable
x_variable_means    = nanmean(x_variables);
%Determine the SD for each x variable
x_variable_SD       = nanstd(x_variables);


%Pull the data for the y variables
y_variables         = numData(:,first_y_column:last_y_column);
%Determine how many y variables there are
number_y_variables  = size(y_variables,2);
%Determine the AVE for each y variable
y_variable_means    = nanmean(y_variables);
%Determine the SD for each y variable
y_variable_SD       = nanstd(y_variables);

%initialize the CorrCoeff (r) and p-value matrices           
r_matrix    = NaN(number_y_variables+1,number_x_variables+1);
pVal_matrix = NaN(number_y_variables+1,number_x_variables+1);


%Loop to correlate each x variable with each y variable
for i = 1:number_x_variables
    xTempVector = x_variables(:,i);          %Pull the ith x variable for the ith round of correlations (redundent...)
    for ii = 1:number_y_variables
        xTempVector = x_variables(:,i);      %Reset xTempVector if any elements were removed in the ii-1th iteration
        yTempVector = y_variables(:,ii);     %Pull the iith y variable for the iith round of correlations
        NaN_index = find(isnan(yTempVector) | isnan(xTempVector));   %Create an index variable that locates NaN in either vector to be correlated
        xTempVector(NaN_index) = [];         %Remove all datapoints in the ith x variable that are NaN in EITHER x or y
        yTempVector(NaN_index) = [];         %Remove all datapoints in the iith y variable that are NaN in EITHER x or y
        
        
        %Calculate and disp the p-value%
        [rho,pVal] = corr(xTempVector,yTempVector);
        

        %Save CorrCoeff (r) to a matrix
        r_matrix(ii+1,i+1)    = rho;
        
        %Save the p-value to a matrix
        pVal_matrix(ii+1,i+1) = pVal;
        
        
    end
    
end

%Convert both of the filled matrices into cell array so that the headers can be added
r_matrix_cell    = num2cell(r_matrix);
pVal_matrix_cell = num2cell(pVal_matrix);


%Add the x variable headers to each of the cell arrays
r_matrix_cell(1,2:number_x_variables+1)    = x_variable_headers;
pVal_matrix_cell(1,2:number_x_variables+1) = x_variable_headers;
%Add the y variable headers to each of the cell arrays
r_matrix_cell(2:number_y_variables+1)    = y_variable_headers;
pVal_matrix_cell(2:number_y_variables+1) = y_variable_headers;


%Create an excel sheet for the CorrCoeffs and p-values (would be nice to put them into the same excel...
writecell(r_matrix_cell,'CorrCoeff_r_matrix.xlsx');
writecell(pVal_matrix_cell,'CorrCoeff_pVal_matrix.xlsx');

fprintf('\n%1.0f "X" variables were correlated with %1.0f "Y" variables\n\n',number_x_variables,number_y_variables)

end












