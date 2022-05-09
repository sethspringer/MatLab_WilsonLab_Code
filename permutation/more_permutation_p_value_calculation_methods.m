
%ONLY NEEDED IF CODE IS STOPPED EARLY
cluster_sum_perm_max_list_shortened = cluster_sum_perm_max_list(1:519);

%Separately per tail

%Positive tail
positive_clusters = cluster_sum_perm_max_list_shortened(cluster_sum_perm_max_list_shortened>0);
positive_clusters_sorted = sort(positive_clusters)


prob_positive_clusters = zeros(length(positive_clusters_sorted),1);

for i = 1:length(positive_clusters_sorted) 

   temp_value = positive_clusters_sorted(i);
   
   prob_positive_clusters(i) = (nnz(positive_clusters_sorted>=temp_value))/length(positive_clusters_sorted);
    

end


plot(positive_clusters_sorted,prob_positive_clusters)


%Now negative
negative_clusters = cluster_sum_perm_max_list_shortened(cluster_sum_perm_max_list_shortened<0);
negative_clusters_sorted = sort(negative_clusters,'descend');


prob_negative_clusters = zeros(length(negative_clusters_sorted),1);

for i = 1:length(negative_clusters_sorted) 

   temp_value = negative_clusters_sorted(i);
   
   prob_negative_clusters(i) = (nnz(negative_clusters_sorted<=temp_value))/length(negative_clusters_sorted);
    

end


plot(negative_clusters_sorted,prob_negative_clusters)













%ABS value method (Bullmore, Sucklin, et al., 1999)
cluster_sum_perm_max_list_shortened_ABS = abs(cluster_sum_perm_max_list_shortened);

sorted_cluster_size_ABS = abs(sorted_cluster_size);

p_values_ABS = zeros(length(sorted_cluster_size),1);

for i = 1:length(sorted_cluster_size_ABS) 

   temp_value = sorted_cluster_size_ABS(i);
   
   p_values_ABS(i) = (nnz(cluster_sum_perm_max_list_shortened_ABS>temp_value))/(length(cluster_sum_perm_max_list_shortened_ABS)+1);
    

end


















%Difference from mean method (Ernst, 2004)

cluster_sum_perm_max_list_shortened_mean = mean(cluster_sum_perm_max_list_shortened);

%Diff scores
cluster_sum_perm_max_list_shortened_mean_diff = abs(cluster_sum_perm_max_list_shortened-cluster_sum_perm_max_list_shortened_mean);
sorted_cluster_size_mean_diff  = abs(sorted_cluster_size-cluster_sum_perm_max_list_shortened_mean);


p_values_mean_diff = zeros(length(sorted_cluster_size),1);

for i = 1:length(sorted_cluster_size_mean_diff) 

   temp_value = sorted_cluster_size_mean_diff(i);
   
   p_values_mean_diff(i) = (nnz(cluster_sum_perm_max_list_shortened_mean_diff>=temp_value))/length(cluster_sum_perm_max_list_shortened_mean_diff);
    

end








