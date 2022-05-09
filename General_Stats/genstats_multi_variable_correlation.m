function genstats_multi_variable_correlation(first_x_column,last_x_column,first_y_column,last_y_column)

%PURPOSE:          Correlate one group of variables with another CAN test for normality and change the test to Spearman accordingly.
%
%REQUIRED INPUTS:  An excel spreadsheet organized with each variable having its own column.
%                  You MUST specify what columns that each variable begins
%                  and ends in...
%
%                  EXAMPLE: If you enter MultipleCorrelations_NormalityTest(3,5,6,8), 
%                  the function will correlate variable 3 with 6, 7, and 8; correlate
%                  variable 4 with 6, 7, and 8; and correlate variable 5 with 6, 7, and 8...
%
%NOTES:            Removes missing values pair-wise.
%                  Requires the "Shapiro-Wilk and Shapiro-Francia normality tests" function.   
%                  This function either tests normality with the S-W or S-F test, depending on the kurtosis of the data...
%                  
%                  Make sure that your data starts in the top left of the excel sheet (A1).
%
%AUTHOR:           Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:  08/31/2020  v1: First working version of program
%                              


[file,path] = uigetfile('*.xlsx','Select the excel you want to work with...');

cd(path)

%Read in the excel file
[numData,txtData,rawData] = xlsread(file);


%Ask the user if they want to run normality testing
prompt       = {'Test for normality? ("Y" or "N")      If no, ONLY Pearson Corr will be used...'};
definput     = {'Y'};
promptAnswer = inputdlg(prompt,'Inputs',[1 35],definput);
%Check to make sure that a valid input was given, otherwise, ask again.
while ~(strcmpi(promptAnswer,'Y') | strcmpi(promptAnswer,'N'))
    promptAnswer = inputdlg(prompt,'Please check your inputs!',[1 50],definput);
end

%Pull the x variable headers
x_variable_headers = rawData(1,first_x_column:last_x_column);
%Pull the y variable headers
y_variable_headers = rawData(1,first_y_column:last_y_column)';

%Pull the data for the x variables
x_variables         = numData(:,first_x_column:last_x_column);
%Determine how many x variables there are
number_x_variables  = size(x_variables,2);


%Pull the data for the y variables
y_variables         = numData(:,first_y_column:last_y_column);
%Determine how many y variables there are
number_y_variables  = size(y_variables,2);

%initialize the CorrCoeff (r) and p-value matrices
r_matrix       = NaN(number_y_variables+1,number_x_variables+1);
pVal_matrix    = NaN(number_y_variables+1,number_x_variables+1);
normality_cell = num2cell(NaN(number_y_variables+1,number_x_variables+1));


%Loop to correlate each x variable with each y variable
for i = 1:number_x_variables
    for ii = 1:number_y_variables
        xTempVector = x_variables(:,i);      %Pull the ith x variable for the ith round of correlations %Reset xTempVector if any elements were removed in the ii-1th iteration
        yTempVector = y_variables(:,ii);     %Pull the iith y variable for the iith round of correlations
        NaN_index = find(isnan(yTempVector) | isnan(xTempVector));   %Create an index variable that locates NaN in either vector to be correlated
        xTempVector(NaN_index) = [];         %Remove all datapoints in the ith x variable that are NaN in EITHER x or y
        yTempVector(NaN_index) = [];         %Remove all datapoints in the iith y variable that are NaN in EITHER x or y
        
        %Checks to see if the user wants to run normality testing
        if strcmpi(promptAnswer,'Y')
            %Runs normality testing for both of the TempVectors
            [H_x,SW_pValue_x,SWstatistic_x] = swtest(xTempVector,0.05);
            [H_y,SW_pValue_y,SWstatistic_y] = swtest(yTempVector,0.05);
            
            
            %If either vector is non-normal, run Spearman. Otherwise, run Pearson
            %Save the type of test used to the cell array (normality_cell)
            if H_x | H_y
                %disp('non-normal')
                [rho,pVal] = corr(xTempVector,yTempVector,'Type','Spearman');
                normality_cell{ii+1,i+1} = 'Non-normal (Spearman)';
            else
                %disp('normal')
                [rho,pVal] = corr(xTempVector,yTempVector);
                normality_cell{ii+1,i+1} = 'Normal (Pearson)';
            end
            
        %If they don't want normality testing, only Pearson Correlation will be run...
        else
            [rho,pVal] = corr(xTempVector,yTempVector);
            normality_cell{ii+1,i+1} = 'Normal (Pearson)';
        end
        
        
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
normality_cell(1,2:number_x_variables+1)   = x_variable_headers;
%Add the y variable headers to each of the cell arrays
r_matrix_cell(2:number_y_variables+1)    = y_variable_headers;
pVal_matrix_cell(2:number_y_variables+1) = y_variable_headers;
normality_cell(2:number_y_variables+1)   = y_variable_headers;

%Write out the excel files for the CorrCoeff matrix, p-values matrix, and Normality of the data/correlation test used.
r_matrix_cell=cell2table(r_matrix_cell(2:end,:),'VariableNames',horzcat('Corr_Coeff',r_matrix_cell(1,2:end)));
writetable(r_matrix_cell,'CorrCoeff_r_values.csv');

pVal_matrix_cell=cell2table(pVal_matrix_cell(2:end,:),'VariableNames',horzcat('p_values',pVal_matrix_cell(1,2:end)));
writetable(pVal_matrix_cell,'CorrCoeff_pValues.csv');

normality_cell=cell2table(normality_cell(2:end,:),'VariableNames',horzcat('Normality_test_and_Correlation_test',normality_cell(1,2:end)));
writetable(normality_cell,'Normality_Testing.csv');

fprintf('\n%1.0f "X" variables were correlated with %1.0f "Y" variables\n\n',number_x_variables,number_y_variables)

end












