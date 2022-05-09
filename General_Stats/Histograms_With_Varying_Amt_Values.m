function Histograms_With_Varying_Amt_Values

%PURPOSE:           Create 5 histograms (with normal distribution curve overlay) of values normally distributed
%                   around zero. Each with more random values.
%                   
%
%REQUIRED INPUTS:   NA
%                   
%		   
%                  
%AUTHOR:            Seth D. Springer, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   08/23/2020  v1: First working version of program.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n100     = randn(1,100);
n1000    = randn(1,1000);
n10000   = randn(1,10000);
n100000  = randn(1,100000);
n1000000 = randn(1,1000000);

figure(1);
histfit(n100);
title('100 Values','FontSize',16,'FontWeight','bold');
movegui(1,[50 550])
xlabel('Std. Deviations From Zero','FontSize',12)
ylabel('Frequency','FontSize',12)

figure(2);
histfit(n1000);
title('1000 Values','FontSize',16,'FontWeight','bold');
movegui(2,[650 550])
xlabel('Std. Deviations From Zero','FontSize',12)
ylabel('Frequency','FontSize',12)

figure(3);
histfit(n10000);
title('10000 Values','FontSize',16,'FontWeight','bold');
movegui(3,[1250 550])
xlabel('Std. Deviations From Zero','FontSize',12)
ylabel('Frequency','FontSize',12)

figure(4);
histfit(n100000);
title('100000 Values','FontSize',16,'FontWeight','bold');
movegui(4,[350 25])
xlabel('Std. Deviations From Zero','FontSize',12)
ylabel('Frequency','FontSize',12)

figure(5);
histfit(n1000000);
title('1000000 Values','FontSize',16,'FontWeight','bold');
movegui(5,[975 25])
xlabel('Std. Deviations From Zero','FontSize',12)
ylabel('Frequency','FontSize',12)


end