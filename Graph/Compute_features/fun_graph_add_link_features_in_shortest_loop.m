function link_feature = fun_graph_add_link_features_in_shortest_loop(link_feature)
% For loop removal classifier 
%
%
%
%

num_link = size(link_feature, 1);
link_feature.shortest_loop_avg_skl_SNR = nan(num_link, 1);
link_feature.shortest_loop_avg_dt_mean = nan(num_link, 1);
link_feature.shortest_loop_avg_int_mean = nan(num_link, 1);
% link_feature.shortest_loop_median_dt_median = nan(num_link, 1);
% link_feature.shortest_loop_median_int_median = nan(num_link, 1);
% Compute loop features
tic
num_link_in_loop_str = numel(tmp_loop_str.link_label);
for iter_link = 1 : num_link_in_loop_str
    if ~isinf(tmp_loop_str.loop_geodesic_length(iter_link)) 
        tmp_current_link_label = tmp_loop_str.link_label(iter_link);
        tmp_loop_link_label = tmp_loop_str.loop_link_label{iter_link};
        assert(~isempty(tmp_loop_link_label), 'There is no link in this loop');
        % Get features 
        tmp_loop_link_lengths = link_feature.length(tmp_loop_link_label);
        tmp_loop_link_weights = tmp_loop_link_lengths ./ tmp_loop_str.loop_length(iter_link);
        tmp_loop_link_num_voxel = link_feature.num_voxel(tmp_loop_link_label);
        tmp_loop_link_skl_SNR = link_feature.skl_SNR(tmp_loop_link_label);
        tmp_loop_link_dt_mean = link_feature.dt_mean(tmp_loop_link_label);
        tmp_loop_link_int_mean = link_feature.int_mean(tmp_loop_link_label);
        
        tmp_not_single_voxel_link_Q = (tmp_loop_link_num_voxel > 1);
        % Exclude single voxel link when its feature is NaN
        link_feature.shortest_loop_avg_skl_SNR(tmp_current_link_label) = tmp_loop_link_skl_SNR(tmp_not_single_voxel_link_Q)' * ...
            tmp_loop_link_weights(tmp_not_single_voxel_link_Q) / (tmp_loop_link_weights' * tmp_not_single_voxel_link_Q);
        
        link_feature.shortest_loop_avg_dt_mean(tmp_current_link_label) = tmp_loop_link_dt_mean' * ...
            tmp_loop_link_weights;
        
        link_feature.shortest_loop_avg_int_mean(tmp_current_link_label) = tmp_loop_link_int_mean' * ...
            tmp_loop_link_weights;    
        
%         tmp_loop_link_dt_median = link_feature.dt_median(tmp_loop_link_label);
%         tmp_loop_link_int_median = link_feature.int_median(tmp_loop_link_label);        
%         link_feature.shortest_loop_median_dt_median(tmp_current_link_label) = median(tmp_loop_link_dt_median, 'omitnan');
%         link_feature.shortest_loop_median_int_median(tmp_current_link_label) = median(tmp_loop_link_int_median, 'omitnan');
    end
end
toc
link_feature.skl_SNR_2_shortest_loop_avg_SNR = link_feature.skl_SNR ./ link_feature.shortest_loop_avg_skl_SNR;
link_feature.dt_mean_2_shortest_loop_avg_dt  = link_feature.dt_mean ./ link_feature.shortest_loop_avg_dt_mean;
link_feature.int_mean_2_shortest_loop_avg_int = link_feature.int_mean ./ link_feature.shortest_loop_avg_int_mean;

end