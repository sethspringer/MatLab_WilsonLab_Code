function seth_permutation_testing_one_way_ANOVA

%This code is demonstrating the example given in the Ernst 2004 paper,
%Table 3, Figures 6 & 7


inputs = {'Number of permutations'};
defaults = {'9999'};
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[n_permutations] = deal(answer{:});
n_permutations = str2num(n_permutations);


n_nii = 14;

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





group1_data = [135;91;111;87;122];
group2_data = [175;130;514;283];
group3_data = [105;147;159;107;194];


data = [group1_data;group2_data;group3_data];
group_labels = {'a';'a';'a';'a';'a';'b';'b';'b';'b';'c';'c';'c';'c';'c'};
group_labels_num = [1,1,1,1,1,2,2,2,2,3,3,3,3,3]';



group1_test = 5*((mean(group1_data))^2)
group2_test = 4*((mean(group2_data))^2)
group3_test = 5*((mean(group3_data))^2)

sum_group_stat = group1_test+group2_test+group3_test;


[p,table,stats,terms] = anovan(data, group_labels_num, 'display','off');
f_stat = str2double(string(table(2,6)));




%Initialize the loading bar
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0);
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);

tic

f_stat_perm_list = zeros(n_permutations,1);
sum_group_stat_perm_list = zeros(n_permutations,1);



%Now I can move on to permuting and sampling clusters...
for permutation_index = 1:n_permutations
    
    
    %I am shuffling the data order
    data_perm = data(valid_perms(permutation_index,:));
    
    
    group1_test_perm = 5*((mean(data_perm(1:5)))^2);
    group2_test_perm = 4*((mean(data_perm(6:9)))^2);
    group3_test_perm = 5*((mean(data_perm(10:end)))^2);
    
    sum_group_stat_perm = group1_test_perm+group2_test_perm+group3_test_perm;
    
    
    [p_perm,table_perm,stats_perm,terms_perm] = anovan(data_perm, group_labels_num, 'display','off');
    f_stat_perm = str2double(string(table_perm(2,6)));
    
    
    f_stat_perm_list(permutation_index) = f_stat_perm;
    sum_group_stat_perm_list(permutation_index) = sum_group_stat_perm;
    
    
    clear f_stat_perm sum_group_stat_perm
    
    if getappdata(wait_bar,'canceling')
        delete(wait_bar)
        error('Permutation testing canceled by user')
    end
    waitbar(permutation_index/n_permutations,wait_bar,sprintf('Permutation #%.0f',permutation_index))
    
end %end of permutation for loop

delete(wait_bar);

total_time = toc





p_val_sum = (nnz(sum_group_stat_perm_list>sum_group_stat))/(n_permutations+1)
p_val_f_stat = (nnz(f_stat_perm_list>f_stat))/(n_permutations+1)










end %end of function