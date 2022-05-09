function genstats_multi_variable_correlation_fisherZ(first_x_column,last_x_column,first_y_column,last_y_column,group1code,group2code)

%PURPOSE:         Correlate multiple variables and compare the resulting correlations between groups (.
%
%REQUIRED INPUTS:  An excel spreadsheet organized with each variable having its 
%                  own column (grouping variable MUST be in the second column).
%
%                  You MUST specify what columns that each variable begins
%                  and ends in...
%                 
%
%                 EXAMPLE: If you enter FisherTransformation_new(3,5,6,8,0,1), 
%                 the function will correlate variable 3 with 6, 7, and 8; correlate
%                 variable 4 with 6, 7, and 8; and correlate variable 5 with 6, 7, and 8... SEPARATELY FOR EACH GROUP
%                 THEN, the function will transform each CorrCoeff (r) into Fisher's Z scores and the Fisher's Z scores will be
%                 compared between groups (answering the question: is the Corr Coeff (r) in question different between your two groups?).
%
%                 It does not matter which group is assigned to "group 1" or "group 2"
%
%
%NOTES:           Removes missing values pair-wise.
%                 The resultant p-values are slightly different than http://vassarstats.net/rdiff.html
%                 bc that website rounds the z to 2 digits while this function rounds to 4 digits.   
%                  
%                 Make sure that your data starts in the top left of the excel sheet (A1).
%
%AUTHOR:           Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:  09/01/2020  v1: First working version of program
%

[file,path] = uigetfile('*.xlsx','Select the excel you want to work with...');

cd(path)

%Read in the excel file
[numData,txtData,rawData] = xlsread(file);


%Pull the x variable headers
x_variable_headers = rawData(1,first_x_column:last_x_column);
%Pull the y variable headers
y_variable_headers = rawData(1,first_y_column:last_y_column)';


%Determine the number of participants in each group
number_in_group_1 = nnz(find(numData(:,2)==group1code));
number_in_group_2 = nnz(find(numData(:,2)==group2code));

row_index_group1 = find(numData(:,2)==group1code);
row_index_group2 = find(numData(:,2)==group2code);

%Pull the data for the x variables
x_variables_group1 = numData(row_index_group1,first_x_column:last_x_column);
x_variables_group2 = numData(row_index_group2,first_x_column:last_x_column);
%Determine how many x variables there are
number_x_variables = size(x_variables_group1,2);


%Pull the data for the y variables
y_variables_group1 = numData(row_index_group1,first_y_column:last_y_column);
y_variables_group2 = numData(row_index_group2,first_y_column:last_y_column);
%Determine how many y variables there are
number_y_variables = size(y_variables_group1,2);



%initialize the CorrCoeff (r), p-value, and element number matrices for each group
r_matrix_group1               = NaN(number_y_variables+1,number_x_variables+1);
r_matrix_group2               = NaN(number_y_variables+1,number_x_variables+1);
pVal_matrix_group1            = NaN(number_y_variables+1,number_x_variables+1);
pVal_matrix_group2            = NaN(number_y_variables+1,number_x_variables+1);
number_elements_matrix_group1 = NaN(number_y_variables+1,number_x_variables+1);
number_elements_matrix_group2 = NaN(number_y_variables+1,number_x_variables+1);


%Loop to correlate each x variable with each y variable within GROUP 1
for i = 1:number_x_variables
    xTempVector = x_variables_group1(:,i);          %Pull the ith x variable for the ith round of correlations (redundent...)
    for ii = 1:number_y_variables
        xTempVector = x_variables_group1(:,i);      %Reset xTempVector if any elements were removed in the ii-1th iteration
        yTempVector = y_variables_group1(:,ii);     %Pull the iith y variable for the iith round of correlations
        NaN_index = find(isnan(yTempVector) | isnan(xTempVector));   %Create an index variable that locates NaN in either vector to be correlated
        xTempVector(NaN_index) = [];         %Remove all datapoints in the ith x variable that are NaN in EITHER x or y
        yTempVector(NaN_index) = [];         %Remove all datapoints in the iith y variable that are NaN in EITHER x or y
        
        n                      = length(xTempVector);
        
        %Calculate and disp the p-value%
        [rho,pVal] = corr(xTempVector,yTempVector);
        

        %Save CorrCoeff (r) to a matrix
        r_matrix_group1(ii+1,i+1)               = rho;
        
        %Save the p-value to a matrix
        pVal_matrix_group1(ii+1,i+1)            = pVal;
        
        %Save the number of elements used in each correlation to a matrix
        number_elements_matrix_group1(ii+1,i+1) = n;
        
    end
    
end

%Loop to correlate each x variable with each y variable within GROUP 2
for i = 1:number_x_variables
    for ii = 1:number_y_variables
        xTempVector            = x_variables_group2(:,i);      %Reset xTempVector if any elements were removed in the ii-1th iteration
        yTempVector            = y_variables_group2(:,ii);     %Pull the iith y variable for the iith round of correlations
        NaN_index              = find(isnan(yTempVector) | isnan(xTempVector));   %Create an index variable that locates NaN in either vector to be correlated
        xTempVector(NaN_index) = [];         %Remove all datapoints in the ith x variable that are NaN in EITHER x or y
        yTempVector(NaN_index) = [];         %Remove all datapoints in the iith y variable that are NaN in EITHER x or y
        
        n                      = length(xTempVector);
        
        %Calculate and disp the p-value%
        [rho,pVal] = corr(xTempVector,yTempVector);
        

        %Save CorrCoeff (r) to a matrix
        r_matrix_group2(ii+1,i+1)               = rho;
        
        %Save the p-value to a matrix
        pVal_matrix_group2(ii+1,i+1)            = pVal;
        
        %Save the number of elements used in each correlation to a matrix
        number_elements_matrix_group2(ii+1,i+1) = n;
        
    end
    
end

%Transform each groups' CorrCoeff (r) matrix into a Fisher's Z matrix
Fishers_z_matrix_group1 = atanh(r_matrix_group1);
Fishers_z_matrix_group2 = atanh(r_matrix_group2);


%Initialize the matrix for the Ztest values (comparing Fisher's Z's of each group
Z_test_matrix = NaN(number_y_variables+1,number_x_variables+1);

for i = 1:number_x_variables
    for ii = 1:number_y_variables
        numerator   = (Fishers_z_matrix_group1(ii+1,i+1)-Fishers_z_matrix_group2(ii+1,i+1));
        denominator = (sqrt((1/(number_elements_matrix_group1(ii+1,i+1)-3))+(1/(number_elements_matrix_group2(ii+1,i+1)-3))));
        Z_test_matrix(ii+1,i+1) = numerator/denominator;
    end
end


% Initialize the matrix that will hold all of the p-values from comparing
% the Z_test values to critical values.
twoTailed_pValues_from_comparing_Z_values = NaN(number_y_variables+1,number_x_variables+1);


for i = 1:number_x_variables
    for ii = 1:number_y_variables
        z = Z_test_matrix(ii+1,i+1);   %Pull out the z_test value in the ii+1th row and the i+1th column.
        
        %The following 10 lines are taken from the "pvaluefromz" function online...
        z = abs(z);
        F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
        p = integral(F, z, 100);
        
        %Commenting this out bc it is unnecessary
        %fprintf ('\nOne tail p-value : %1.6f', p);
        %fprintf ('\nTwo tail p-value : %1.6f\n' ,p*2)
        
        %Return two tail p-value
        pvalue = p*2;
        
        twoTailed_pValues_from_comparing_Z_values(ii+1,i+1) = pvalue;
        
    end
end

%Convert the filled matrices into cell array so that the headers can be added
r_matrix_cell_group1            = num2cell(r_matrix_group1);
r_matrix_cell_group2            = num2cell(r_matrix_group2);
pVal_matrix_cell_group1         = num2cell(pVal_matrix_group1);
pVal_matrix_cell_group2         = num2cell(pVal_matrix_group2);
number_elements_cell_group1     = num2cell(number_elements_matrix_group1);
number_elements_cell_group2     = num2cell(number_elements_matrix_group2);

FisherZ_Comparison_pValues_cell = num2cell(twoTailed_pValues_from_comparing_Z_values);

%Add the x variable headers to each of the cell arrays
r_matrix_cell_group1(1,2:number_x_variables+1)            = x_variable_headers;
r_matrix_cell_group2(1,2:number_x_variables+1)            = x_variable_headers;
pVal_matrix_cell_group1(1,2:number_x_variables+1)         = x_variable_headers;
pVal_matrix_cell_group2(1,2:number_x_variables+1)         = x_variable_headers;
number_elements_cell_group1(1,2:number_x_variables+1)     = x_variable_headers;
number_elements_cell_group2(1,2:number_x_variables+1)     = x_variable_headers;

FisherZ_Comparison_pValues_cell(1,2:number_x_variables+1) = x_variable_headers;

%Add the y variable headers to each of the cell arrays
r_matrix_cell_group1(2:number_y_variables+1)            = y_variable_headers;
r_matrix_cell_group2(2:number_y_variables+1)            = y_variable_headers;
pVal_matrix_cell_group1(2:number_y_variables+1)         = y_variable_headers;
pVal_matrix_cell_group2(2:number_y_variables+1)         = y_variable_headers;
number_elements_cell_group1(2:number_y_variables+1)     = y_variable_headers;
number_elements_cell_group2(2:number_y_variables+1)     = y_variable_headers;

FisherZ_Comparison_pValues_cell(2:number_y_variables+1) = y_variable_headers;


%Convert the cells into tables and write them out as excel (.csv) files
r_matrix_table_group1=cell2table(r_matrix_cell_group1(2:end,:),'VariableNames',horzcat('Corr_Coeff_group1',r_matrix_cell_group1(1,2:end)));
writetable(r_matrix_table_group1,'CorrCoeff_r_values_group1.csv');
r_matrix_table_group2=cell2table(r_matrix_cell_group2(2:end,:),'VariableNames',horzcat('Corr_Coeff_group2',r_matrix_cell_group2(1,2:end)));
writetable(r_matrix_table_group2,'CorrCoeff_r_values_group2.csv');

%Convert the cells into tables and write out the pValues a/w each correlation ran
pVal_matrix_table_group1=cell2table(pVal_matrix_cell_group1(2:end,:),'VariableNames',horzcat('CorrCoeff_p_values_group1',pVal_matrix_cell_group1(1,2:end)));
writetable(pVal_matrix_table_group1,'CorrCoeff_pValues_group1.csv');
pVal_matrix_table_group2=cell2table(pVal_matrix_cell_group2(2:end,:),'VariableNames',horzcat('CorrCoeff_p_values_group2',pVal_matrix_cell_group2(1,2:end)));
writetable(pVal_matrix_table_group2,'CorrCoeff_pValues_group2.csv');

%Convert the cells into tables and write out the number of elements in each correlation
number_elements_table_group1=cell2table(number_elements_cell_group1(2:end,:),'VariableNames',horzcat('Number_of_Subjects_group1',number_elements_cell_group1(1,2:end)));
writetable(number_elements_table_group1,'Numb_Subjects_n_group1.csv');
number_elements_table_group2=cell2table(number_elements_cell_group2(2:end,:),'VariableNames',horzcat('Number_of_Subjects_group2',number_elements_cell_group2(1,2:end)));
writetable(number_elements_table_group2,'Numb_Subjects_n_group2.csv');

%Convert the cells into tables and write out the pValues a/w the Fisher's Z comparisons
FisherZ_Comparison_pValues_table=cell2table(FisherZ_Comparison_pValues_cell(2:end,:),'VariableNames',horzcat('FisherZ_GroupComparison_pValues',FisherZ_Comparison_pValues_cell(1,2:end)));
writetable(FisherZ_Comparison_pValues_table,'FisherZ_GroupComparison_pValues.csv');

fprintf('\n%1.0f "X" variables were correlated with %1.0f "Y" variables and these correlations were compared between groups!\n\n',number_x_variables,number_y_variables)

end
