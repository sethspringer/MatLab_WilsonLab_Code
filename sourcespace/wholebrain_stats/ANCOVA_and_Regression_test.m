
y_easy = [12,32,12,54,65,34,45,87]';

y_hard = [23,43,65,34,65,68,87,56]';

group = [-1,-1,-1,-1,-1,1,1,1]';

cov   = [32,14,64,2,23,64,7,3]';



t_ANOVA = table(group,y_easy,y_hard,'VariableNames',{'Group','Easy','Hard'});


rm_table_ANOVA = table([1 2]','VariableNames',{'Condition'});

rm_ANOVA = fitrm(t_ANOVA,'Hard-Easy~Group','WithinDesign',rm_table_ANOVA);

model_ANOVA = ranova(rm_ANOVA);



%With covariate
t_ANCOVA = table(group,y_easy,y_hard,cov,...
'VariableNames',{'Group','Easy','Hard','Covariate'});
rm_table_ANCOVA = table([1 2]','VariableNames',{'Condition'});

rm_ANCOVA = fitrm(t_ANCOVA,'Hard-Easy~Group+Covariate','WithinDesign',rm_table_ANCOVA);

model_ANCOVA = ranova(rm_ANCOVA);





%Group effect
y_both = [y_easy,y_hard];
y_ave = mean(y_both,2);

tbl = table(group,y_ave,...
'VariableNames',{'Group','Condition_Ave'});


lm = fitlm(tbl,'Condition_Ave~Group')


%Group effect, with covariate
y_both = [y_easy,y_hard];
y_ave = mean(y_both,2);

tbl_cov = table(group,y_ave,cov,...
'VariableNames',{'Group','Condition_Ave','Covariate'});


lm_cov = fitlm(tbl_cov,'Condition_Ave~Group+Covariate')








