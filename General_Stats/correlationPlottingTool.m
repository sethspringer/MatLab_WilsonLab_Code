function correlationPlottingTool(x_column,y_column)

%PURPOSE:          Read in one x and one y variable from a spreadsheet and plot them together in a scatterplot. 
%                  Display best fit line.
%
%REQUIRED INPUTS:  Variable columns within the excel.
%                  
%
%AUTHOR:           Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:  09/03/2020  v1: First working version of program
%                              


[file,path] = uigetfile('*.xlsx','Select the excel you want to work with...');

cd(path)

%Read in the excel file
[numData,txtData,rawData] = xlsread(file);


%Pull the x variable headers
x_variable_header = rawData(1,x_column);
%Pull the y variable headers
y_variable_header = rawData(1,y_column);


%Pull the data for the x variables
x_variable = numData(:,x_column);

%Pull the data for the y variables
y_variable = numData(:,y_column);


NaN_index = find(isnan(y_variable) | isnan(x_variable));   %Create an index variable that locates NaN in either vector to be correlated
x_variable(NaN_index) = [];                               %Remove all datapoints in the x variable that are NaN in EITHER x or y
y_variable(NaN_index) = [];                               %Remove all datapoints in the y variable that are NaN in EITHER x or y


%Produce scatterplot with bounds determined by the input vectors
scatter(x_variable,y_variable)
xlim([(min(x_variable)-(max(x_variable)*0.1)) (max(x_variable)+(max(x_variable)*0.1))])
ylim([(min(y_variable)-(max(y_variable)*0.1)) (max(y_variable)+(max(y_variable)*0.1))])
lsline                          %Adds a least-squares (ls; best-fit) line to the graph.

xlabel(x_variable_header)
ylabel(y_variable_header)

end












