function timeseries_montecarlo_toolbox_seth_test
%PURPOSE:           Perform cluster-based permutation sampling of
%                   timeseries data
%
%REQUIRED INPUTS:   (1) MonteCarloPermutation_ToolboxWorksheet with correct
%                   parameters specified and data
%                   (2) Data specifications, dependent on test
%		        
%AUTHOR:            Alex I. Wiesman, DICoN Lab, University of Nebraska Medical Center
%VERSION HISTORY:   11/29/2017  v1: First working version of program
%                   10/08/2019  v2: (1) Now uses maximum cluster from
%                   each permutation (more robust to Type I error - also more stringent) 
%                   (2) compiled all functions into one *.m file (3)
%                   uses only non-overlapping permutations (4) now
%                   only runs </= maximum number of possible permutations
%                   (5) no longer clusters opposing effects that are
%                   continuous
%                   11/13/2019 v3: (1) Fixed bug for unpaired
%                   t-tests, where unequal group sizes would previously
%                   error-out

%prompt for data file and read initial parameters
[data_file,path,~] = uigetfile('*.xlsx','Please Select Worksheet','MonteCarloPermutation_ToolboxWorksheet');
cd(path)
[~,test_type,~] = xlsread(data_file,'Parameters','B2');
sig_cutoff = xlsread(data_file,'Parameters','B5');
num_permutations = xlsread(data_file,'Parameters','B8');
tails = xlsread(data_file,'Parameters','B11');

if isempty(test_type) | isempty(sig_cutoff) | isempty(num_permutations)
    error('Not all test parameters have been specified properly. Refer to notes on the Worksheet.')
end

if strcmp(test_type,'correlation') 
    timeseries_montecarlo_toolbox_correlation(data_file,sig_cutoff,num_permutations,tails)
elseif strcmp(test_type,'baseline-t')
    timeseries_montecarlo_toolbox_baselinet(data_file,sig_cutoff,num_permutations,tails)
elseif strcmp(test_type,'paired-t')
    timeseries_montecarlo_toolbox_pairedt(data_file,sig_cutoff,num_permutations,tails)
elseif strcmp(test_type,'unpaired-t')
    timeseries_montecarlo_toolbox_unpairedt(data_file,sig_cutoff,num_permutations,tails)
elseif strcmp(test_type,'anova')
    timeseries_montecarlo_toolbox_oneway_anova(data_file,sig_cutoff,num_permutations,tails)
else
    error('No valid statistical test specified.')
end

function  timeseries_montecarlo_toolbox_baselinet(data_file,sig_cutoff,num_permutations,tails)

%collect additional test-specific data%
inputs = {'Epoch Start (ms)','Time Sampling Resolution (ms)','Baseline Start (ms)','Baseline End (ms)'};
defaults = {'0','25','0','0'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[epoch_start,time_samp,base_start,base_end] = deal(answer{:});
time_samp =str2num(time_samp);
epoch_start =str2num(epoch_start);
base_start =str2num(base_start);
base_end =str2num(base_end);

%read time series data from spreadsheet%
timeseries_data = xlsread(data_file,'Timeseries Data');

%calculate time stamps based on user input%
total_time = time_samp*(size(timeseries_data,2)-1);
base_start_index = ((base_start-epoch_start)/time_samp)+1;
base_end_index = ((base_end-epoch_start)/time_samp)+1;
end_time = total_time+epoch_start;
time = epoch_start:time_samp:end_time;

%randomly shuffle the baseline data across the length of the time series%
baseline_repetitions = ceil(size(timeseries_data,2)/(base_end_index-base_start_index));
for i = 1:size(timeseries_data,1)
    base = repmat(timeseries_data(i,base_start_index:base_end_index),1,baseline_repetitions);
    rand_base = base(randperm(size(base,2)));
    indiv_baselines(i,:) = rand_base;
end
indiv_baselines = indiv_baselines(:,1:size(timeseries_data,2));

%Compute max number of valid permutations (within a reasonable time limit)%
num_ts = size(timeseries_data,1);
for i = 1:num_permutations
    perm_counter = 1;
    while 1
        if perm_counter > num_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',num_permutations,max_perms);
            num_permutations = size(valid_perms,1);
            break
        end
        new_perm = randi([0 1],[1,num_ts]);
        perm_counter = perm_counter+1;
        if exist('valid_perms','var') && ~isempty(valid_perms)
            if ~ismember(new_perm,valid_perms,'rows') & sum(new_perm) ~= 0 & sum(new_perm) ~= num_ts
                valid_perms(i,:) = new_perm;
                    break
            end
        elseif sum(new_perm) ~= 0 & sum(new_perm) ~= num_ts
                valid_perms(i,:) = new_perm;
                break
        end
    end
end

valid_perms = logical(valid_perms)';

%calculate the initial timeseries of coefficients
[~,sig_timeseries,~,stats] = ttest(timeseries_data,indiv_baselines);
coeff_timeseries = stats.tstat;
if tails == 1
    sig_timeseries = sig_timeseries/2;
end

%find all values that exceed significance cutoff and determine which are contiguous with other significant values
candidate_clusters = find(sig_timeseries < sig_cutoff/100);
consec_bins = diff(candidate_clusters) == 1;
if size(candidate_clusters,2) > 0
    if size(candidate_clusters,2) == 1
       consec_bins(1) = 0;
    elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
        consec_bins(end+1) = 0;
    end
else
    error('No potentially significant clusters identified prior to permutation testing. Consider raising the significance level.')
end

%grow each cluster of significant values to determine sampling cutoffs
bin_num = 1;
bin_start = [];
bin_end = [];
bin_ref = 1;
while 1
    if bin_ref > length(candidate_clusters)
       break
    end
    if bin_ref > length(consec_bins)
        break
    end
    if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
        if consec_bins(1,bin_ref) 
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                while 1
                    if bin_ref+1 <= length(consec_bins) 
                        if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                        elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                            break
                        else
                            break
                        end
                    else 
                        break
                    end
                end
        else
            bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
            bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
        end
            bin_ref = find(candidate_clusters == bin_end(1,bin_num));
            bin_num = bin_num+1;
    else
        bin_ref = bin_ref+1;
    end
end

%compute sum of test coefficients for each cluster
cluster_coefficient_sum = zeros(1,length(bin_start));
for i = 1:length(bin_start)
    cluster_coefficient_sum(1,i) = sum(coeff_timeseries(1,bin_start(1,i):bin_end(1,i)));
end

cluster_numbers = (1:size(bin_start,2))';
cluster_times(:,1) = time(bin_start);
cluster_times(:,2) = time(bin_end);
final_coeff_timeseries = coeff_timeseries;
average_timeseries = mean(timeseries_data,1);

%for each permutation, shuffle the data across participants then compute a new timeseries of coefficients
perm_coefficient_sum = zeros(num_permutations,1);
coeff_timeseries = zeros(1,size(timeseries_data,2));
sig_timeseries = zeros(1,size(timeseries_data,2));
full_data = [timeseries_data;indiv_baselines];
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0)
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
for i = 1:num_permutations
    randgroup1 = timeseries_data;
    randgroup1(valid_perms(:,i),:) = indiv_baselines(valid_perms(:,i),:);
    randgroup2 = indiv_baselines;
    randgroup2(valid_perms(:,i),:) = timeseries_data(valid_perms(:,i),:);

    [~,sig_timeseries,~,stats] = ttest(randgroup1,randgroup2);
    coeff_timeseries = stats.tstat;

if tails == 1
    sig_timeseries = sig_timeseries/2;
end
    
%identify temporally contiguous clusters of significance
    candidate_clusters = find(sig_timeseries < sig_cutoff/100);
    consec_bins = diff(candidate_clusters) == 1;
    if size(candidate_clusters,2) > 0
        if size(candidate_clusters,2) == 1
           consec_bins(1) = 0;
        elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
            consec_bins(end+1) = 0;
        end
    end
        bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end
    
%sum the coefficients within each cluster, and build a distribution of these sums
    if ~isempty(bin_start)
        bin_start(bin_start == 0) = [];
        bin_end(bin_end == 0) = [];
        for ii = 1:length(bin_start)
            bin_sums(ii) = sum(coeff_timeseries(1,bin_start(1,ii):bin_end(1,ii)));
        end
        [~,max_bin_sum_index] = max(abs(bin_sums));
        perm_coefficient_sum(i) = bin_sums(max_bin_sum_index);
        clear bin_sums;
    end

if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Computation Canceled')
end
waitbar(i/num_permutations,wait_bar,sprintf('Permutation #%.0f',i))
end
delete(wait_bar);

positive_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum >= 0));
negative_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum <= 0));

%calculate the p-value of each cluster from the original timeseries
cluster_pvals = zeros(1,size(cluster_coefficient_sum,2));
for i = 1:size(cluster_coefficient_sum,2) 
    if cluster_coefficient_sum(1,i)>0
        cluster_pvals(1,i) = 1-((sum(sum(positive_perm_coefficient_sum<cluster_coefficient_sum(1,i)))/sum(sum(~isnan(positive_perm_coefficient_sum)))));
    else
        cluster_pvals(1,i) = 1-((sum(sum(negative_perm_coefficient_sum>cluster_coefficient_sum(1,i)))/sum(sum(~isnan(negative_perm_coefficient_sum)))));
    end
end

if tails == 1
    cluster_pvals = cluster_pvals/2;
end

%write the results to a new worksheet in the original spreadsheet
final_data = [cluster_numbers, cluster_times, cluster_pvals'];
xlswrite(data_file,[{'Cluster #'},{'Start (ms)'},{'End (ms)'},{'p-value'},{'# Permutations'}],'Results');
xlswrite(data_file,num_permutations,'Results','E2');
xlswrite(data_file,final_data,'Results','A2');

%plot the average timeseries and r-value timeseries data, and highlight clusters (red = statistically significant at original threshold, blue = non-significant)
subplot(2,1,1)
plot(time,final_coeff_timeseries);
ylabel('Test Statistic (t)');
subplot(2,1,2)
plot(time,average_timeseries);
ylabel('Average Timeseries');
hold on
xlabel('Time (ms)');
for i = 1:size(final_data,1)
    if cluster_pvals(1,i) < (sig_cutoff/100)
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min(average_timeseries) min(average_timeseries) max(average_timeseries) max(average_timeseries)],'r','EdgeColor','r','FaceAlpha',.3)
    else
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min(average_timeseries) min(average_timeseries) max(average_timeseries) max(average_timeseries)],'b','EdgeColor','b','FaceAlpha', .3)  
    end
end

%plot the distribution of the permuted null hypothesis, and highlight the significant clusters
sorted_perm_clusters = sort(perm_coefficient_sum);
unique_sorted_clusters = unique(sort(perm_coefficient_sum));
hist_ranges = linspace(unique_sorted_clusters(1,1),unique_sorted_clusters(end,1),(num_permutations/10));
hist_counts = zeros(1,size(hist_ranges,2));
for i = 1:size(hist_ranges,2)-1
    hist_counts(1,i) = sum(sorted_perm_clusters<hist_ranges(1,i+1) & sorted_perm_clusters>=hist_ranges(1,i));       
end
figure
bar(hist_ranges,hist_counts);
xlabel('Cluster Sum');
ylabel('# of Observations');
hold on

for i = 1:size(cluster_coefficient_sum,2)
    if cluster_pvals(1,i)<(sig_cutoff/100)
        line([cluster_coefficient_sum(1,i) cluster_coefficient_sum(1,i)],[0 max(hist_counts)],'Color','red','LineStyle','--');
    end
end
        
function  timeseries_montecarlo_toolbox_correlation(data_file,sig_cutoff,num_permutations,tails)

inputs = {'Epoch Start (ms)','Time Sampling Resolution (ms)','# of Covariates (of No Interest)'};
defaults = {'0','25','0'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[epoch_start,time_samp,num_covars] = deal(answer{:});
time_samp =str2num(time_samp);
epoch_start =str2num(epoch_start);
num_covars =str2num(num_covars);

timeseries_data = xlsread(data_file,'Timeseries Data');
orig_timeseries_data = timeseries_data;
covar_interest = xlsread(data_file,'Covariate of Interest');
static_covars_nointerest = xlsread(data_file,'Covariate(s) of No Interest');
num_time_covars = num_covars - size(static_covars_nointerest,2);

if num_time_covars == 0 & ~num_covars == 0
    covar_type = 'static';
    if ~isequal(size(static_covars_nointerest,1),size(timeseries_data,1))
        error('Covariates of no interest must have the same dimensions as the timeseries!')
    end
elseif num_time_covars == num_covars & ~num_covars == 0
    covar_type = 'time';
elseif ~num_covars == 0
    covar_type = 'both';
    if ~isequal(size(static_covars_nointerest,1),size(timeseries_data,1))
        error('Covariates of no interest must have the same dimensions as the timeseries!')
    end
else 
    covar_type = [];
end

if num_time_covars == 1
        covars_nointerest1 = xlsread(data_file,'Timeseries Data 2');
        if ~isequal(size(covars_nointerest1),size(timeseries_data))
            error('Covariates of no interest must have the same dimensions as the timeseries!')
        end
        clear covars_nointerest
        covars_nointerest = permute(covars_nointerest1,[1 3 2]);
elseif num_time_covars == 2
        covars_nointerest1 = xlsread(data_file,'Timeseries Data 2');
        covars_nointerest2 = xlsread(data_file,'Timeseries Data 3');
        if ~isequal(size(covars_nointerest1),size(covars_nointerest2),size(timeseries_data))
            error('Covariates of no interest must have the same dimensions as the timeseries!')
        end
        clear covars_nointerest
        covars_nointerest = cat(3,covars_nointerest1,covars_nointerest2);
        covars_nointerest = permute(covars_nointerest,[1 3 2]);
elseif num_time_covars == 3
        covars_nointerest1 = xlsread(data_file,'Timeseries Data 2');
        covars_nointerest2 = xlsread(data_file,'Timeseries Data 3');
        covars_nointerest3 = xlsread(data_file,'Timeseries Data 4');
        if ~isequal(size(covars_nointerest1),size(covars_nointerest2),size(covars_nointerest3),size(timeseries_data))
            error('Covariates of no interest must have the same dimensions as the timeseries!')
        end
        clear covars_nointerest
        covars_nointerest = cat(3,covars_nointerest1,covars_nointerest2,covars_nointerest3);
        covars_nointerest = permute(covars_nointerest,[1 3 2]);
elseif num_time_covars > 3
    error('No more than three temporally-extended (timeseries) covariates of no interest are accepted at this time!')
end

%Compute max number of valid permutations (within a reasonable time limit)
num_ts = size(timeseries_data,1);
for i = 1:num_permutations
    perm_counter = 1;
    while 1
        if perm_counter > num_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',num_permutations,max_perms);
            num_permutations = size(valid_perms,1);
            break
        end
        new_perm = randperm(num_ts);   %Create a random permutation ordering
        perm_counter = perm_counter+1; 
        if exist('valid_perms','var') && ~isempty(valid_perms) %If valid_perms is a variable that is not empty
            if ~ismember(new_perm,valid_perms,'rows') & ~isequal(new_perm,[1:num_ts]) %Before adding in the new random permutation vector, make sure that it differs from existing ones and is not linearly arranged
                valid_perms(i,:) = new_perm;
                break
            end
        elseif ~isequal(new_perm,[1:num_ts]) %Make sure that the new random order is not linearly arranged
            valid_perms(i,:) = new_perm;
        end
    end
end

total_time = time_samp*(size(timeseries_data,2)-1);
end_time = total_time+epoch_start;
time = epoch_start:time_samp:end_time;

%calculate the initial timeseries of coefficients
coeff_timeseries = zeros(1,size(timeseries_data,2));
sig_timeseries = zeros(1,size(timeseries_data,2));
    if num_covars == 0
        for i = 1:size(timeseries_data,2)
            [r,p] = corrcoef(timeseries_data(:,i),covar_interest);
            coeff_timeseries(1,i) = r(1,2);
            sig_timeseries(1,i) = p(1,2);
        end
    elseif strcmp(covar_type,'static')
        for i = 1:size(timeseries_data,2)
            [coeff_timeseries(1,i),sig_timeseries(1,i)] = partialcorr(timeseries_data(:,i),covar_interest,static_covars_nointerest);
        end
    elseif strcmp(covar_type,'time')
        for i = 1:size(timeseries_data,2)
            [coeff_timeseries(1,i),sig_timeseries(1,i)] = partialcorr(timeseries_data(:,i),covar_interest,covars_nointerest(:,:,i));        
        end
    elseif strcmp(covar_type,'both')
        for i = 1:size(timeseries_data,2)
            [coeff_timeseries(1,i),sig_timeseries(1,i)] = partialcorr(timeseries_data(:,i),covar_interest,[covars_nointerest(:,:,i),static_covars_nointerest]);        
        end  
    end
    
if tails == 1
    sig_timeseries = sig_timeseries/2;
end

%find all values that exceed significance cutoff and determine which are contiguous with other significant values
candidate_clusters = find(sig_timeseries < sig_cutoff/100);
consec_bins = diff(candidate_clusters) == 1;
if size(candidate_clusters,2) > 0
    if size(candidate_clusters,2) == 1
       consec_bins(1) = 0;
    elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
        consec_bins(end+1) = 0;
    end
else
    error('No potentially significant clusters identified prior to permutation testing. Consider raising the significance level.')
end

%grow each cluster of significant values to determine sampling cutoffs
bin_num = 1;
bin_start = [];
bin_end = [];
bin_ref = 1;
while 1
    if bin_ref > length(candidate_clusters)
       break
    end
    if bin_ref > length(consec_bins)
        break
    end
    if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
        if consec_bins(1,bin_ref) 
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                while 1
                    if bin_ref+1 <= length(consec_bins) 
                        if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                        elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                            break
                        else
                            break
                        end
                    else 
                        break
                    end
                end
        else
            bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
            bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
        end
            bin_ref = find(candidate_clusters == bin_end(1,bin_num));
            bin_num = bin_num+1;
    else
        bin_ref = bin_ref+1;
    end
end


%compute sum of test coefficients for each cluster
cluster_coefficient_sum = zeros(1,length(bin_start));
for i = 1:length(bin_start)
    cluster_coefficient_sum(1,i) = sum(coeff_timeseries(1,bin_start(1,i):bin_end(1,i)));
end

cluster_numbers = (1:size(bin_start,2))';
cluster_times(:,1) = time(bin_start);
cluster_times(:,2) = time(bin_end);
final_coeff_timeseries = coeff_timeseries;
average_timeseries = mean(timeseries_data,1);

%for each permutation, shuffle the data across participants then compute a new timeseries of coefficients
perm_coefficient_sum = zeros(num_permutations,1);
coeff_timeseries = zeros(1,size(timeseries_data,2));
sig_timeseries = zeros(1,size(timeseries_data,2));
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0);
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
for i = 1:num_permutations
    timeseries_data = orig_timeseries_data(valid_perms(i,:),:); %Randomly shuffling the timeseries data across participants
    
    if num_covars == 0
        for ii = 1:size(timeseries_data,2)
            [r,p] = corrcoef(timeseries_data(:,ii),covar_interest);
            coeff_timeseries(1,ii) = r(1,2);
            sig_timeseries(1,ii) = p(1,2);
        end
    elseif strcmp(covar_type,'static')
        for ii = 1:size(timeseries_data,2)
            [coeff_timeseries(1,ii),sig_timeseries(1,ii)] = partialcorr(timeseries_data(:,ii),covar_interest,static_covars_nointerest);
        end
    elseif strcmp(covar_type,'time')
        for ii = 1:size(timeseries_data,2)
            [coeff_timeseries(1,ii),sig_timeseries(1,ii)] = partialcorr(timeseries_data(:,ii),covar_interest,covars_nointerest(:,:,ii));           
        end
    elseif strcmp(covar_type,'both')
        for ii = 1:size(timeseries_data,2)
            [coeff_timeseries(1,ii),sig_timeseries(1,ii)] = partialcorr(timeseries_data(:,ii),covar_interest,[covars_nointerest(:,:,ii),static_covars_nointerest]);        
        end              
    end
  
if tails == 1
    sig_timeseries = sig_timeseries/2;
end

%identify temporally contiguous clusters of significance
    candidate_clusters = find(sig_timeseries < sig_cutoff/100);
    consec_bins = diff(candidate_clusters) == 1;
    if size(candidate_clusters,2) > 0
        if size(candidate_clusters,2) == 1
           consec_bins(1) = 0;
        elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
            consec_bins(end+1) = 0;
        end
    end
    bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end
    
%sum the coefficients within each cluster, and build a distribution of these sums
    if ~isempty(bin_start)
        bin_start(bin_start == 0) = [];
        bin_end(bin_end == 0) = [];
        for ii = 1:length(bin_start)
            bin_sums(ii) = sum(coeff_timeseries(1,bin_start(1,ii):bin_end(1,ii)));
        end
        [~,max_bin_sum_index] = max(abs(bin_sums));
        perm_coefficient_sum(i) = bin_sums(max_bin_sum_index); %Max cluster values for the i th permutation
        clear bin_sums;
    end

if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Computation Canceled')
end
waitbar(i/num_permutations,wait_bar,sprintf('Permutation #%.0f',i))
end
delete(wait_bar);

positive_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum >= 0)); %grab positive clusters
negative_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum <= 0)); %grab negative clusters

%calculate the p-value of each cluster from the original timeseries %%%SETH i am here
cluster_pvals = zeros(1,size(cluster_coefficient_sum,2));
for i = 1:size(cluster_coefficient_sum,2) 
    if cluster_coefficient_sum(1,i)>0
        cluster_pvals(1,i) = 1-((sum(sum(positive_perm_coefficient_sum<cluster_coefficient_sum(1,i)))/sum(sum(~isnan(positive_perm_coefficient_sum)))));
    else
        cluster_pvals(1,i) = 1-((sum(sum(negative_perm_coefficient_sum>cluster_coefficient_sum(1,i)))/sum(sum(~isnan(negative_perm_coefficient_sum)))));
    end
end

if tails == 1
    cluster_pvals = cluster_pvals/2;
end

%write the results to a new worksheet in the original spreadsheet
final_data = [cluster_numbers, cluster_times, cluster_pvals'];
xlswrite(data_file,[{'Cluster #'},{'Start (ms)'},{'End (ms)'},{'p-value'},{'# Permutations'}],'Results');
xlswrite(data_file,num_permutations,'Results','E2');
xlswrite(data_file,final_data,'Results','A2');

%plot the average timeseries and r-value timeseries data, and highlight clusters (red = statistically significant at original threshold, blue = non-significant)
subplot(2,1,1)
plot(time,final_coeff_timeseries);
ylabel('Test Statistic (r)');
subplot(2,1,2)
plot(time,average_timeseries);
ylabel('Average Timeseries');
hold on
xlabel('Time (ms)');
for i = 1:size(final_data,1)
    if cluster_pvals(1,i) < (sig_cutoff/100)
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min(average_timeseries) min(average_timeseries) max(average_timeseries) max(average_timeseries)],'r','EdgeColor','r','FaceAlpha',.3)
    else
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min(average_timeseries) min(average_timeseries) max(average_timeseries) max(average_timeseries)],'b','EdgeColor','b','FaceAlpha', .3)  
    end
end

%plot the distribution of the permuted null hypothesis, and highlight the significant clusters
sorted_perm_clusters = sort(perm_coefficient_sum);
unique_sorted_clusters = unique(sort(perm_coefficient_sum));
hist_ranges = linspace(unique_sorted_clusters(1,1),unique_sorted_clusters(end,1),(num_permutations/10));
hist_counts = zeros(1,size(hist_ranges,2));
for i = 1:size(hist_ranges,2)-1
    hist_counts(1,i) = sum(sorted_perm_clusters<hist_ranges(1,i+1) & sorted_perm_clusters>=hist_ranges(1,i));       
end
figure
bar(hist_ranges,hist_counts);
xlabel('Cluster Sum');
ylabel('# of Observations');
hold on

for i = 1:size(cluster_coefficient_sum,2)
    if cluster_pvals(1,i)<(sig_cutoff/100)
        line([cluster_coefficient_sum(1,i) cluster_coefficient_sum(1,i)],[0 max(hist_counts)],'Color','red','LineStyle','--');
    end
end

function  timeseries_montecarlo_toolbox_pairedt(data_file,sig_cutoff,num_permutations,tails)

inputs = {'Epoch Start (ms)','Time Sampling Resolution (ms)'};
defaults = {'0','25'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[epoch_start,time_samp] = deal(answer{:});
time_samp =str2num(time_samp);
epoch_start =str2num(epoch_start);

timeseries_data1 = xlsread(data_file,'Timeseries Data');
timeseries_data2 = xlsread(data_file,'Timeseries Data 2');

if ~size(timeseries_data1) == size(timeseries_data2)
    error('Both sets of timeseries data must have the same dimensions when computing a paired-samples test!')
end

%Compute max number of valid permutations (within a reasonable time limit)
num_ts = size(timeseries_data1,1);
for i = 1:num_permutations
    perm_counter = 1;
    while 1
        if perm_counter > num_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',num_permutations,max_perms);
            num_permutations = size(valid_perms,1);
            break
        end
        new_perm = randi([0 1],[1,num_ts]);
        perm_counter = perm_counter+1;
        if exist('valid_perms','var') && ~isempty(valid_perms)
            if ~ismember(new_perm,valid_perms,'rows') & sum(new_perm) ~= 0 & sum(new_perm) ~= num_ts
                valid_perms(i,:) = new_perm;
                    break
            end
        elseif sum(new_perm) ~= 0 & sum(new_perm) ~= num_ts
                valid_perms(i,:) = new_perm;
                break
        end
    end
end

valid_perms = logical(valid_perms)';

total_time = time_samp*(size(timeseries_data1,2)-1);
end_time = total_time+epoch_start;
time = epoch_start:time_samp:end_time;

%calculate the initial timeseries of coefficients
[~,sig_timeseries,~,stats] = ttest(timeseries_data1,timeseries_data2);
coeff_timeseries = stats.tstat;

if tails == 1
    sig_timeseries = sig_timeseries/2;
end

%find all values that exceed significance cutoff and determine which are contiguous with other significant values
candidate_clusters = find(sig_timeseries < sig_cutoff/100);
consec_bins = diff(candidate_clusters) == 1;
if size(candidate_clusters,2) > 0
    if size(candidate_clusters,2) == 1
       consec_bins(1) = 0;
    elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
        consec_bins(end+1) = 0;
    end
else
    error('No potentially significant clusters identified prior to permutation testing. Consider raising the significance level.')
end

%grow each cluster of significant values to determine sampling cutoffs
    bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end

%compute sum of test coefficients for each cluster
cluster_coefficient_sum = zeros(1,length(bin_start));
for i = 1:length(bin_start)
    cluster_coefficient_sum(1,i) = sum(coeff_timeseries(1,bin_start(1,i):bin_end(1,i)));
end

cluster_numbers = (1:size(bin_start,2))';
cluster_times(:,1) = time(bin_start);
cluster_times(:,2) = time(bin_end);
final_coeff_timeseries = coeff_timeseries;
average_timeseries1 = mean(timeseries_data1,1);
average_timeseries2 = mean(timeseries_data2,1);

%for each permutation, shuffle the data across participants then compute a new timeseries of coefficients
perm_coefficient_sum = zeros(num_permutations,1);
coeff_timeseries = zeros(1,size(timeseries_data1,2));
sig_timeseries = zeros(1,size(timeseries_data1,2));
% full_data = [timeseries_data1;timeseries_data2];
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0)
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
for i = 1:num_permutations
    randgroup1 = timeseries_data1;
    randgroup1(valid_perms(:,i),:) = timeseries_data2(valid_perms(:,i),:);
    randgroup2 = timeseries_data2;
    randgroup2(valid_perms(:,i),:) = timeseries_data1(valid_perms(:,i),:);

    [~,sig_timeseries,~,stats] = ttest(randgroup1,randgroup2);
    coeff_timeseries = stats.tstat;
    
if tails == 1
    sig_timeseries = sig_timeseries/2;
end
    
%identify temporally contiguous clusters of significance
    candidate_clusters = find(sig_timeseries < sig_cutoff/100);
    consec_bins = diff(candidate_clusters) == 1;
    if size(candidate_clusters,2) > 0
        if size(candidate_clusters,2) == 1
           consec_bins(1) = 0;
        elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
            consec_bins(end+1) = 0;
        end
    end
            bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end
    
%sum the coefficients within each cluster, and build a distribution of these sums
    if ~isempty(bin_start)
        bin_start(bin_start == 0) = [];
        bin_end(bin_end == 0) = [];
        for ii = 1:length(bin_start)
            bin_sums(ii) = sum(coeff_timeseries(1,bin_start(1,ii):bin_end(1,ii)));
        end
        [~,max_bin_sum_index] = max(abs(bin_sums));
        perm_coefficient_sum(i) = bin_sums(max_bin_sum_index);
        clear bin_sums;
    end

if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Computation Canceled')
end
waitbar(i/num_permutations,wait_bar,sprintf('Permutation #%.0f',i))
end
delete(wait_bar);

positive_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum >= 0));
negative_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum <= 0));

%calculate the p-value of each cluster from the original timeseries
cluster_pvals = zeros(1,size(cluster_coefficient_sum,2));
for i = 1:size(cluster_coefficient_sum,2) 
    if cluster_coefficient_sum(1,i)>0
        cluster_pvals(1,i) = 1-((sum(sum(positive_perm_coefficient_sum<cluster_coefficient_sum(1,i)))/sum(sum(~isnan(positive_perm_coefficient_sum)))));
    else
        cluster_pvals(1,i) = 1-((sum(sum(negative_perm_coefficient_sum>cluster_coefficient_sum(1,i)))/sum(sum(~isnan(negative_perm_coefficient_sum)))));
    end
end

if tails == 1
    cluster_pvals = cluster_pvals/2;
end

%write the results to a new worksheet in the original spreadsheet
final_data = [cluster_numbers, cluster_times, cluster_pvals'];
xlswrite(data_file,[{'Cluster #'},{'Start (ms)'},{'End (ms)'},{'p-value'},{'# Permutations'}],'Results');
xlswrite(data_file,num_permutations,'Results','E2');
xlswrite(data_file,final_data,'Results','A2');

%plot the average timeseries and r-value timeseries data, and highlight clusters (red = statistically significant at original threshold, blue = non-significant)
subplot(2,1,1)
plot(time,final_coeff_timeseries);
ylabel('Test Statistic (t)');
subplot(2,1,2)
plot(time,average_timeseries1,time,average_timeseries2);
ylabel('Average Timeseries');
hold on
xlabel('Time (ms)');
for i = 1:size(final_data,1)
    if cluster_pvals(1,i) < (sig_cutoff/100)
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2]) min([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2])],'r','EdgeColor','r','FaceAlpha',.3)
    else
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2]) min([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2])],'b','EdgeColor','b','FaceAlpha', .3)  
    end
end

%plot the distribution of the permuted null hypothesis, and highlight the significant clusters
sorted_perm_clusters = sort(perm_coefficient_sum);
unique_sorted_clusters = unique(sort(perm_coefficient_sum));
hist_ranges = linspace(unique_sorted_clusters(1,1),unique_sorted_clusters(end,1),(num_permutations/10));
hist_counts = zeros(1,size(hist_ranges,2));
for i = 1:size(hist_ranges,2)-1
    hist_counts(1,i) = sum(sorted_perm_clusters<hist_ranges(1,i+1) & sorted_perm_clusters>=hist_ranges(1,i));       
end
figure
bar(hist_ranges,hist_counts);
xlabel('Cluster Sum');
ylabel('# of Observations');
hold on

for i = 1:size(cluster_coefficient_sum,2)
    if cluster_pvals(1,i)<(sig_cutoff/100)
        line([cluster_coefficient_sum(1,i) cluster_coefficient_sum(1,i)],[0 max(hist_counts)],'Color','red','LineStyle','--');
    end
end

function  timeseries_montecarlo_toolbox_unpairedt(data_file,sig_cutoff,num_permutations,tails)

inputs = {'Epoch Start (ms)','Time Sampling Resolution (ms)'};
defaults = {'0','25'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[epoch_start,time_samp] = deal(answer{:});
time_samp =str2num(time_samp);
epoch_start =str2num(epoch_start);

timeseries_data1 = xlsread(data_file,'Timeseries Data');
timeseries_data2 = xlsread(data_file,'Timeseries Data 2');

if ~size(timeseries_data1,2) == size(timeseries_data2,2)
    error('Both sets of timeseries data must have the same number of samples in the time domain!')
end

%Compute max number of valid permutations (within a reasonable time limit)
initial_assignments = zeros(1,(size(timeseries_data1,1)+size(timeseries_data2,1)));
initial_assignments(1,1:size(timeseries_data1,1)) = 1;
initial_assignments(1,size(timeseries_data1,1)+1:size(timeseries_data1,1)+size(timeseries_data2,1)) = 2;
num_ts = size(timeseries_data1,1)+size(timeseries_data2,1);
for i = 1:num_permutations
    perm_counter = 1;
    while 1
        if perm_counter > num_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',num_permutations,max_perms);
            num_permutations = size(valid_perms,1);
            break
        end
        new_perm = initial_assignments(randperm(size(initial_assignments,2)));
        perm_counter = perm_counter+1;
        if exist('valid_perms','var') && ~isempty(valid_perms)
            if ~ismember(new_perm,valid_perms,'rows') & ~((mean(new_perm(1,1:size(timeseries_data1,1))) == 1) && (mean(new_perm(1,size(timeseries_data1,1)+1:size(timeseries_data2,1))) == 2))
                valid_perms(i,:) = new_perm;
                    break
            end
        elseif ~((mean(new_perm(1,1:size(timeseries_data1,1))) == 1) && (mean(new_perm(1,size(timeseries_data1,1)+1:size(timeseries_data2,1))) == 2))
                valid_perms(i,:) = new_perm;
                break
        end
    end
end

total_time = time_samp*(size(timeseries_data1,2)-1);
end_time = total_time+epoch_start;
time = epoch_start:time_samp:end_time;

%calculate the initial timeseries of coefficients
[~,sig_timeseries,~,stats] = ttest2(timeseries_data1,timeseries_data2);
coeff_timeseries = stats.tstat;

if tails == 1
    sig_timeseries = sig_timeseries/2;
end

%find all values that exceed significance cutoff and determine which are contiguous with other significant values
candidate_clusters = find(sig_timeseries < sig_cutoff/100);
consec_bins = diff(candidate_clusters) == 1;
if size(candidate_clusters,2) > 0
    if size(candidate_clusters,2) == 1
       consec_bins(1) = 0;
    elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
        consec_bins(end+1) = 0;
    end
else
    error('No potentially significant clusters identified prior to permutation testing. Consider raising the significance level.')
end

%grow each cluster of significant values to determine sampling cutoffs
bin_num = 1;
bin_start = [];
bin_end = [];
bin_ref = 1;
while 1
    if bin_ref > length(candidate_clusters)
       break
    end
    if bin_ref > length(consec_bins)
        break
    end
    if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
        if consec_bins(1,bin_ref) 
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                while 1
                    if bin_ref+1 <= length(consec_bins) 
                        if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                        elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                            break
                        else
                            break
                        end
                    else 
                        break
                    end
                end
        else
            bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
            bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
        end
            bin_ref = find(candidate_clusters == bin_end(1,bin_num));
            bin_num = bin_num+1;
    else
        bin_ref = bin_ref+1;
    end
end

if length(candidate_clusters) == 1
    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
    bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
end

%compute sum of test coefficients for each cluster
cluster_coefficient_sum = zeros(1,length(bin_start));
for i = 1:length(bin_start)
    cluster_coefficient_sum(1,i) = sum(coeff_timeseries(1,bin_start(1,i):bin_end(1,i)));
end

cluster_numbers = (1:size(bin_start,2))';
cluster_times(:,1) = time(bin_start);
cluster_times(:,2) = time(bin_end);
final_coeff_timeseries = coeff_timeseries;
average_timeseries1 = mean(timeseries_data1,1);
average_timeseries2 = mean(timeseries_data2,1);

%for each permutation, shuffle the data across participants then compute a new timeseries of coefficients
perm_coefficient_sum = zeros(num_permutations,1);
coeff_timeseries = zeros(1,size(timeseries_data1,2));
sig_timeseries = zeros(1,size(timeseries_data1,2));
full_data = [timeseries_data1;timeseries_data2];
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0)
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
for i = 1:num_permutations
    randgroup1 = full_data(find(valid_perms(i,:)==1),:);
    randgroup2 = full_data(find(valid_perms(i,:)==2),:);
%     randgroup1 = timeseries_data1;
%     randgroup1(valid_perms(:,i),:) = timeseries_data2(valid_perms(:,i),:);
%     randgroup2 = timeseries_data2;
%     randgroup2(valid_perms(:,i),:) = timeseries_data1(valid_perms(:,i),:);

    [~,sig_timeseries,~,stats] = ttest2(randgroup1,randgroup2);
    coeff_timeseries = stats.tstat;
    
    
if tails == 1
    sig_timeseries = sig_timeseries/2;
end
    
%identify temporally contiguous clusters of significance
    candidate_clusters = find(sig_timeseries < sig_cutoff/100);
    consec_bins = diff(candidate_clusters) == 1;
    if size(candidate_clusters,2) > 0
        if size(candidate_clusters,2) == 1
           consec_bins(1) = 0;
        elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
            consec_bins(end+1) = 0;
        end
    end
    bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end
    
%sum the coefficients within each cluster, and build a distribution of these sums
    if ~isempty(bin_start)
        bin_start(bin_start == 0) = [];
        bin_end(bin_end == 0) = [];
        for ii = 1:length(bin_start)
            bin_sums(ii) = sum(coeff_timeseries(1,bin_start(1,ii):bin_end(1,ii)));
        end
        [~,max_bin_sum_index] = max(abs(bin_sums));
        perm_coefficient_sum(i) = bin_sums(max_bin_sum_index);
        clear bin_sums;
    end

if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Computation Canceled')
end
waitbar(i/num_permutations,wait_bar,sprintf('Permutation #%.0f',i))
end
delete(wait_bar);

positive_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum > 0));
negative_perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum < 0));

%calculate the p-value of each cluster from the original timeseries
cluster_pvals = zeros(1,size(cluster_coefficient_sum,2));
for i = 1:size(cluster_coefficient_sum,2) 
    if cluster_coefficient_sum(1,i)>0
        cluster_pvals(1,i) = 1-((sum(sum(positive_perm_coefficient_sum<cluster_coefficient_sum(1,i)))/sum(sum(~isnan(positive_perm_coefficient_sum)))));
    else
        cluster_pvals(1,i) = 1-((sum(sum(negative_perm_coefficient_sum>cluster_coefficient_sum(1,i)))/sum(sum(~isnan(negative_perm_coefficient_sum)))));
    end
end

if tails == 1
    cluster_pvals = cluster_pvals/2;
end
save

%write the results to a new worksheet in the original spreadsheet
final_data = [cluster_numbers, cluster_times, cluster_pvals'];
xlswrite(data_file,[{'Cluster #'},{'Start (ms)'},{'End (ms)'},{'p-value'},{'# Permutations'}],'Results');
xlswrite(data_file,num_permutations,'Results','E2');
xlswrite(data_file,final_data,'Results','A2');

%plot the average timeseries and r-value timeseries data, and highlight clusters (red = statistically significant at original threshold, blue = non-significant)
subplot(2,1,1)
plot(time,final_coeff_timeseries);
ylabel('Test Statistic (t)');
subplot(2,1,2)
plot(time,average_timeseries1,time,average_timeseries2);
ylabel('Average Timeseries');
hold on
xlabel('Time (ms)');
for i = 1:size(final_data,1)
    if cluster_pvals(1,i) < (sig_cutoff/100)
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2]) min([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2])],'r','EdgeColor','r','FaceAlpha',.3)
    else
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2]) min([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2]) max([average_timeseries1 average_timeseries2])],'b','EdgeColor','b','FaceAlpha', .3)  
    end
end

%plot the distribution of the permuted null hypothesis, and highlight the significant clusters
sorted_perm_clusters = sort([negative_perm_coefficient_sum; positive_perm_coefficient_sum]);
hist_ranges = linspace(sorted_perm_clusters(1,1),sorted_perm_clusters(end,1),(num_permutations/10));
hist_counts = zeros(1,size(hist_ranges,2));
for i = 1:size(hist_ranges,2)-1
    hist_counts(1,i) = sum(sorted_perm_clusters<hist_ranges(1,i+1) & sorted_perm_clusters>=hist_ranges(1,i));       
end
figure
bar(hist_ranges,hist_counts);
xlabel('Cluster Sum');
ylabel('# of Observations');
hold on

for i = 1:size(cluster_coefficient_sum,2)
    if cluster_pvals(1,i)<(sig_cutoff/100)
        line([cluster_coefficient_sum(1,i) cluster_coefficient_sum(1,i)],[0 max(hist_counts)],'Color','red','LineStyle','--');
    end
end

function  timeseries_montecarlo_toolbox_oneway_anova(data_file,sig_cutoff,num_permutations,tails)

inputs = {'Epoch Start (ms)','Time Sampling Resolution (ms)'};
defaults = {'0','25'};	
answer = inputdlg(inputs, 'Please Input Parameters', 2, defaults,'on');
[epoch_start,time_samp] = deal(answer{:});
time_samp =str2num(time_samp);
epoch_start =str2num(epoch_start);

timeseries_data1 = xlsread(data_file,'Timeseries Data');
timeseries_data2 = xlsread(data_file,'Timeseries Data 2');
timeseries_data3 = xlsread(data_file,'Timeseries Data 3');

if ~size(timeseries_data1,2) == size(timeseries_data2,2) == size(timeseries_data3,2)
    error('All sets of timeseries data must have the same number of samples in the time domain!')
end

initial_assignments = zeros(1,(size(timeseries_data1,1)+size(timeseries_data2,1)+size(timeseries_data3,1)));
initial_assignments(1,1:size(timeseries_data1,1)) = 1;
initial_assignments(1,size(timeseries_data1,1)+1:size(timeseries_data1,1)+size(timeseries_data2,1)) = 2;
initial_assignments(1,size(timeseries_data1,1)+size(timeseries_data2,1)+1:size(timeseries_data1,1)+size(timeseries_data2,1)+size(timeseries_data3,1)) = 3;

%Compute max number of valid permutations (within a reasonable time limit)
num_ts = size(timeseries_data1,1)+size(timeseries_data2,1)+size(timeseries_data3,1);
for i = 1:num_permutations
    perm_counter = 1;
    while 1
        if perm_counter > num_permutations
            warning('Could not compute %d unique permutations given these data - running with maximum of %d iterations!',num_permutations,max_perms);
            num_permutations = size(valid_perms,1);
            break
        end
        new_perm = initial_assignments(randperm(size(initial_assignments,2)));
        perm_counter = perm_counter+1;
        if exist('valid_perms','var') && ~isempty(valid_perms)
            if ~ismember(new_perm,valid_perms,'rows') & ~((mean(new_perm(1,1:size(timeseries_data1,1))) == 1) && (mean(new_perm(1,size(timeseries_data1,1)+1:size(timeseries_data2,1))) == 2) && (mean(new_perm(1,size(timeseries_data2,1)+1:size(timeseries_data3,1))) == 3)) 
                valid_perms(i,:) = new_perm;
                    break
            end
        elseif ~((mean(new_perm(1,1:size(timeseries_data1,1))) == 1) && (mean(new_perm(1,size(timeseries_data1,1)+1:size(timeseries_data2,1))) == 2) && (mean(new_perm(1,size(timeseries_data2,1)+1:size(timeseries_data3,1))) == 3)) 
                valid_perms(i,:) = new_perm;
                break
        end
    end
end

groups = cell(1,(size(timeseries_data1,1)+size(timeseries_data2,1)+size(timeseries_data3,1)));
groups(1,1:size(timeseries_data1,1)) = {'A'};
groups(1,size(timeseries_data1,1)+1:size(timeseries_data1,1)+size(timeseries_data2,1)) = {'B'};
groups(1,size(timeseries_data1,1)+size(timeseries_data2,1)+1:size(timeseries_data1,1)+size(timeseries_data2,1)+size(timeseries_data3,1)) = {'C'};

total_time = time_samp*(size(timeseries_data1,2)-1);
end_time = total_time+epoch_start;
time = epoch_start:time_samp:end_time;

%calculate the initial timeseries of coefficients
sig_timeseries = zeros(1,size(timeseries_data1,2));
coeff_timeseries = zeros(1,size(timeseries_data1,2));
for i = 1:size(timeseries_data1,2)
[sig_timeseries(1,i),tbl,~] = anova1([timeseries_data1(:,i); timeseries_data2(:,i); timeseries_data3(:,i)]', groups,'off');
coeff_timeseries(1,i) = tbl{2,5};
end

%find all values that exceed significance cutoff and determine which are contiguous with other significant values
candidate_clusters = find(sig_timeseries < sig_cutoff/100);
consec_bins = diff(candidate_clusters) == 1;
if size(candidate_clusters,2) > 0
    if size(candidate_clusters,2) == 1
       consec_bins(1) = 0;
    elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
        consec_bins(end+1) = 0;
    end
else
    error('No potentially significant clusters identified prior to permutation testing. Consider raising the significance level.')
end

%grow each cluster of significant values to determine sampling cutoffs
bin_num = 1;
bin_start = [];
bin_end = [];
bin_ref = 1;
while 1
    if bin_ref > length(candidate_clusters)
       break
    end
    if bin_ref > length(consec_bins)
        break
    end
    if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
        if consec_bins(1,bin_ref) 
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                while 1
                    if bin_ref+1 <= length(consec_bins) 
                        if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                        elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                            bin_ref = bin_ref+1;
                            bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                            break
                        else
                            break
                        end
                    else 
                        break
                    end
                end
        else
            bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
            bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
        end
            bin_ref = find(candidate_clusters == bin_end(1,bin_num));
            bin_num = bin_num+1;
    else
        bin_ref = bin_ref+1;
    end
end

%compute sum of test coefficients for each cluster
cluster_coefficient_sum = zeros(1,length(bin_start));
for i = 1:length(bin_start)
    cluster_coefficient_sum(1,i) = sum(coeff_timeseries(1,bin_start(1,i):bin_end(1,i)));
end

cluster_numbers = (1:size(bin_start,2))';
cluster_times(:,1) = time(bin_start);
cluster_times(:,2) = time(bin_end);
final_coeff_timeseries = coeff_timeseries;
average_timeseries1 = mean(timeseries_data1,1);
average_timeseries2 = mean(timeseries_data2,1);
average_timeseries3 = mean(timeseries_data3,1);

%for each permutation, shuffle the data across participants then compute a new timeseries of coefficients
perm_coefficient_sum = zeros(num_permutations,1);
coeff_timeseries = zeros(1,size(timeseries_data1,2));
sig_timeseries = zeros(1,size(timeseries_data1,2));
full_data = [timeseries_data1;timeseries_data2;timeseries_data3];
wait_bar = waitbar(0,'1','Name','Running test...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait_bar,'canceling',0)
wait_bar_obj=findobj(wait_bar,'Type','Patch');
set(wait_bar_obj,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
for i = 1:num_permutations
    randgroup1 = full_data(find(valid_perms(i,:)==1),:);
    randgroup2 = full_data(find(valid_perms(i,:)==2),:);
    randgroup3 = full_data(find(valid_perms(i,:)==3),:);
       
    for ii = 1:size(timeseries_data1,2)
        [sig_timeseries(1,ii),tbl,~] = anova1([randgroup1(:,ii); randgroup2(:,ii); randgroup3(:,ii)]', groups,'off');
        coeff_timeseries(1,ii) = tbl{2,5};
    end
    
%identify temporally contiguous clusters of significance
    candidate_clusters = find(sig_timeseries < sig_cutoff/100);
    consec_bins = diff(candidate_clusters) == 1;
    if size(candidate_clusters,2) > 0
        if size(candidate_clusters,2) == 1
           consec_bins(1) = 0;
        elseif (candidate_clusters(1,end)-candidate_clusters(1,end-1)) ~= 1
            consec_bins(end+1) = 0;
        end
    end
        bin_num = 1;
    bin_start = [];
    bin_end = [];
    bin_ref = 1;
    while 1
        if bin_ref > length(candidate_clusters)
           break
        end
        if bin_ref > length(consec_bins)
            break
        end
        if ~ismember(candidate_clusters(1,bin_ref),bin_start) && ~ismember(candidate_clusters(1,bin_ref),bin_end)
            if consec_bins(1,bin_ref) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                    bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                    bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                    while 1
                        if bin_ref+1 <= length(consec_bins) 
                            if consec_bins(1,bin_ref+1) & isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref+1);
                            elseif consec_bins(1,bin_ref+1) & ~isequal(sign(coeff_timeseries(candidate_clusters(1,bin_ref))),sign(coeff_timeseries(candidate_clusters(1,bin_ref+1))))
                                bin_ref = bin_ref+1;
                                bin_end(1,bin_num) = candidate_clusters(1,bin_ref-1);
                                break
                            else
                                break
                            end
                        else 
                            break
                        end
                    end
            else
                bin_start(1,bin_num) = candidate_clusters(1,bin_ref);
                bin_end(1,bin_num) = candidate_clusters(1,bin_ref);
            end
                bin_ref = find(candidate_clusters == bin_end(1,bin_num));
                bin_num = bin_num+1;
        else
            bin_ref = bin_ref+1;
        end
    end
    
%sum the coefficients within each cluster, and build a distribution of these sums
    if ~isempty(bin_start)
        bin_start(bin_start == 0) = [];
        bin_end(bin_end == 0) = [];
        for ii = 1:length(bin_start)
            bin_sums(ii) = sum(coeff_timeseries(1,bin_start(1,ii):bin_end(1,ii)));
        end
        [~,max_bin_sum_index] = max(abs(bin_sums));
        perm_coefficient_sum(i) = bin_sums(max_bin_sum_index);
        clear bin_sums;
    end

if getappdata(wait_bar,'canceling')
    delete(wait_bar)
    error('Computation Canceled')
end
waitbar(i/num_permutations,wait_bar,sprintf('Permutation #%.0f',i))
end
delete(wait_bar);

perm_coefficient_sum = (perm_coefficient_sum(perm_coefficient_sum >= 0));

%calculate the p-value of each cluster from the original timeseries
cluster_pvals = zeros(1,size(cluster_coefficient_sum,2));
for i = 1:size(cluster_coefficient_sum,2) 
        cluster_pvals(1,i) = 1-((sum(sum(perm_coefficient_sum<cluster_coefficient_sum(1,i)))/sum(sum(~isnan(perm_coefficient_sum)))));
end

if tails == 1
    warning('Two-tailed tests cannot be run on ANOVA summary statistics! Running one-tailed instead.')
end

%write the results to a new worksheet in the original spreadsheet
final_data = [cluster_numbers, cluster_times, cluster_pvals'];
xlswrite(data_file,[{'Cluster #'},{'Start (ms)'},{'End (ms)'},{'p-value'},{'# Permutations'}],'Results');
xlswrite(data_file,num_permutations,'Results','E2');
xlswrite(data_file,final_data,'Results','A2');

%plot the average timeseries and r-value timeseries data, and highlight clusters (red = statistically significant at original threshold, blue = non-significant)
subplot(2,1,1)
plot(time,final_coeff_timeseries);
ylabel('Test Statistic (F)');
subplot(2,1,2)
plot(time,average_timeseries1,time,average_timeseries2,time,average_timeseries3);
ylabel('Average Timeseries');
hold on
xlabel('Time (ms)');
for i = 1:size(final_data,1)
    if cluster_pvals(1,i) < (sig_cutoff/100)
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2 average_timeseries3]) min([average_timeseries1 average_timeseries2 average_timeseries3]) max([average_timeseries1 average_timeseries2 average_timeseries3]) max([average_timeseries1 average_timeseries2 average_timeseries3])],'r','EdgeColor','r','FaceAlpha',.3)
    else
    patch([cluster_times(i,1) cluster_times(i,2) cluster_times(i,2) cluster_times(i,1)], [min([average_timeseries1 average_timeseries2 average_timeseries3]) min([average_timeseries1 average_timeseries2 average_timeseries3]) max([average_timeseries1 average_timeseries2 average_timeseries3]) max([average_timeseries1 average_timeseries2 average_timeseries3])],'b','EdgeColor','b','FaceAlpha', .3)  
    end
end

%plot the distribution of the permuted null hypothesis, and highlight the significant clusters
sorted_perm_clusters = sort(perm_coefficient_sum);
unique_sorted_clusters = unique(sort(perm_coefficient_sum));
hist_ranges = linspace(unique_sorted_clusters(1,1),unique_sorted_clusters(end,1),(num_permutations/10));
hist_counts = zeros(1,size(hist_ranges,2));
for i = 1:size(hist_ranges,2)-1
    hist_counts(1,i) = sum(sorted_perm_clusters<hist_ranges(1,i+1) & sorted_perm_clusters>=hist_ranges(1,i));       
end
figure
bar(hist_ranges,hist_counts);
xlabel('Cluster Sum');
ylabel('# of Observations');
hold on

for i = 1:size(cluster_coefficient_sum,2)
    if cluster_pvals(1,i)<(sig_cutoff/100)
        line([cluster_coefficient_sum(1,i) cluster_coefficient_sum(1,i)],[0 max(hist_counts)],'Color','red','LineStyle','--');
    end
end