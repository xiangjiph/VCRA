function new_link_feature = fun_graph_get_link_features_in_shortest_loop_by_link_label(link_feature, link_label)
% For loop removal classifier 
%
%
%
%
num_link = size(link_label, 1);
new_link_feature = link_feature(link_label, :);
new_link_feature.shortest_loop_avg_skl_SNR = nan(num_link, 1);
new_link_feature.shortest_loop_avg_dt_mean = nan(num_link, 1);
new_link_feature.shortest_loop_avg_int_mean = nan(num_link, 1);
new_link_feature.shortest_loop_median_dt_median = nan(num_link, 1);
new_link_feature.shortest_loop_median_int_median = nan(num_link, 1);
% Compute loop features
for iter_link = 1 : num_link
    if ~isinf(new_link_feature.shortest_loop_geodesic_lenght(iter_link)) 
        tmp_loop_link_label = new_link_feature.shortest_loop_link_label{iter_link};
        assert(~isempty(tmp_loop_link_label), 'There is no link in this loop');
        % Get features 
        tmp_loop_link_lengths = link_feature.length(tmp_loop_link_label);
        tmp_loop_link_weights = tmp_loop_link_lengths ./ new_link_feature.shortest_loop_lenght(iter_link);
        tmp_loop_link_num_voxel = link_feature.num_voxel(tmp_loop_link_label);
        tmp_loop_link_skl_SNR = link_feature.skl_SNR(tmp_loop_link_label);
        tmp_loop_link_dt_mean = link_feature.dt_mean(tmp_loop_link_label);
        tmp_loop_link_int_mean = link_feature.int_mean(tmp_loop_link_label);
        
        tmp_not_single_voxel_link_Q = (tmp_loop_link_num_voxel > 1);
        % Exclude single voxel link when its feature is NaN
        new_link_feature.shortest_loop_avg_skl_SNR(iter_link) = tmp_loop_link_skl_SNR(tmp_not_single_voxel_link_Q)' * ...
            tmp_loop_link_weights(tmp_not_single_voxel_link_Q) / (tmp_loop_link_weights' * tmp_not_single_voxel_link_Q);
        
        new_link_feature.shortest_loop_avg_dt_mean(iter_link) = tmp_loop_link_dt_mean' * ...
            tmp_loop_link_weights;
        
        new_link_feature.shortest_loop_avg_int_mean(iter_link) = tmp_loop_link_int_mean' * ...
            tmp_loop_link_weights;    
        
        tmp_loop_link_dt_median = link_feature.dt_median(tmp_loop_link_label);
        tmp_loop_link_int_median = link_feature.int_median(tmp_loop_link_label);        
        new_link_feature.shortest_loop_median_dt_median(iter_link) = median(tmp_loop_link_dt_median, 'omitnan');
        new_link_feature.shortest_loop_median_int_median(iter_link) = median(tmp_loop_link_int_median, 'omitnan');
    end
end
new_link_feature.skl_SNR_2_shortest_loop_avg_SNR = new_link_feature.skl_SNR ./ new_link_feature.shortest_loop_avg_skl_SNR;
new_link_feature.dt_mean_2_shortest_loop_avg_dt  = new_link_feature.dt_mean ./ new_link_feature.shortest_loop_avg_dt_mean;
new_link_feature.int_mean_2_shortest_loop_avg_int = new_link_feature.int_mean ./ new_link_feature.shortest_loop_avg_int_mean;
end