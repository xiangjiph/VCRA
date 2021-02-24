function result_str = fun_analysis_region_length_radius_distributions(link_features_T)

% Length distribution
segment_data = struct;
segment_data.radius = link_features_T.dt_median;
segment_data.length = link_features_T.length;
is_valid_Q = isfinite(segment_data.radius) & isfinite(segment_data.length);
if ~all(is_valid_Q)
    fprintf('%d out of %d links with nonfinite radius or length have been removed\n', ...
        nnz(~is_valid_Q), numel(is_valid_Q));
end
segment_data.length = segment_data.length(is_valid_Q);
segment_data.radius = segment_data.radius(is_valid_Q);
% Derivative features
segment_data.surface_area = 2 * pi * segment_data.length .* segment_data.radius;
segment_data.volume = pi * segment_data.length .* segment_data.radius .^ 2;

field_name = fieldnames(segment_data);
% Basic statistics 
result_str = struct;
for iter_features = 1 : numel(field_name)
    tmp_fn = field_name{iter_features};
    tmp_data = segment_data.(tmp_fn);
    result_str.dist.(tmp_fn) = fun_analysis_get_basic_statistics(tmp_data);
end
%% Compute CLF directly from CDF
length_sorted = sort(segment_data.length, 'ascend');
tmp_cum_length = cumsum(length_sorted);
[tmp_unique_length, tmp_idx, ~] = unique(length_sorted, 'stable');
tmp_idx = [tmp_idx(2:end) - 1; numel(length_sorted)];
tmp_cum_length = tmp_cum_length(tmp_idx);
cum_cap_length_itp = griddedInterpolant(tmp_unique_length, tmp_cum_length);
result_str.cumulate_length.x = tmp_unique_length;
result_str.cumulate_length.y_itp = cum_cap_length_itp;
end