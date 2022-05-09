function CorrelationCalculationAndPlotV2(x,y)

%PURPOSE:           Calculate the correlation coeff (r) for the two input vectors (x and y)
%                   and create a scatterplot of x and y.
%
%                   This function uses the equation II.A on page 15 of the correlation book... 
%
%REQUIRED INPUTS:   The two input vectors (x and y) must be held within the
%                   workspace. The inputs MUST be ROW vectors.
%		   
%
%
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   08/23/2020  v1: First working version of program.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Converting the vector inputs into column vectors (the "fit" f(x) requires column vectors as inputs)...
x = reshape(x,[],1);
y = reshape(y,[],1);


n                    = size(x,1);            %Determine the size of the input vectors.
df                   = n-1;                  %Calculate the degrees of freedom.
sumArrayNumerator    = zeros(n,1);           %Preallocate a vector that with hold the numerator values that need to be summed.
sumArrayDenominator1 = zeros(n,1);           %Preallocate a vector that with hold the first set of denominator values that need to be summed.
sumArrayDenominator2 = zeros(n,1);           %Preallocate a vector that with hold the second set of denominator values that need to be summed.

%Descriptive statistics%
meanY   = mean(y);
meanX   = mean(x);
stdDevY = std(y);
stdDevX = std(x);

%Calculate z-scores for each array element and store them in the array "sumArray"%
for i = 1:n
    sumArrayNumerator(i)    = (x(i)-meanX)*(y(i)-meanY);
    sumArrayDenominator1(i) = ((x(i)-meanX).^2);
    sumArrayDenominator2(i) = ((y(i)-meanY).^2);
end

disp(' ')                       %Blank line to separate the function input from the outputs in the Command Window.


%Calculate and disp the r value%
r = (sum(sumArrayNumerator))/(sqrt(((sum(sumArrayDenominator1))*(sum(sumArrayDenominator2)))));
CorrCoeffText = ['CorrCoeff (r)              = ',num2str(r)];
disp(CorrCoeffText)

%Calculate and disp the r^2 value%
CoeffOfDetermination = r^2;
CorrCoeffSquaredText = ['CoeffOfDetermination (r^2) = ',num2str(CoeffOfDetermination)];
disp(CorrCoeffSquaredText)

%Calculate and disp the t-stat%
tStat = (r*(sqrt((n-2))))/(sqrt((1-(r^2))));
tStatText = ['t-statistic                = ',num2str(tStat)];
disp(tStatText)

%Calculate and disp the p-value%
[rho,pVal] = corr(x,y);
pValText = ['p-value                    = ',num2str(pVal)];
disp(pValText)


disp(' ')                      %Blank line the separate the slope of the line from the other values.

%Calculate and disp the slope of the best-fit line%
bestFitLine = fit(x,y,'poly1');                  %The fit function determines the best-fit equation. poly1 = linear fit.
SlopeText = ['Slope of the best-fit line is ',num2str(bestFitLine.p1),'!'];
disp(SlopeText)


disp(' ')                      %Blank line the separate the end of the outputs from the next input line.


%Produce scatterplot with bounds determined by the input vectors
scatter(x,y)
xlim([(min(x)-(max(x)*0.1)) (max(x)+(max(x)*0.1))])
ylim([(min(y)-(max(y)*0.1)) (max(y)+(max(y)*0.1))])
lsline                          %Adds a least-squares (ls; best-fit) line to the graph (same as regression line...)

axis square

%Add lines to mark the means of X and Y, respectively%
xline(meanX);
yline(meanY);

end