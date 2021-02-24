function registered_str = fun_match_point_cloud_by_cpd(pc_fix_sub, pc_mov_sub, opt)

[registered_str, ~] = cpd_register(pc_fix_sub, pc_mov_sub, opt);
dist_1_to_2t = pdist2(pc_fix_sub, registered_str.Y);

switch opt.establish_correspondance_method
    case {'min_dist', 'minimum_distance'}
        [selected_dist, moving_matched_idx] = min(dist_1_to_2t, [], 2);
        if ~isfield(opt, 'max_pair_dist')
            opt.max_pair_dist = 4;
            warning('Maximum pair distance is not specificed. Default value: %f\n', ...
                opt.max_pair_dist);
        end
        fixed_matched_idx = find(selected_dist <= opt.max_pair_dist);
        moving_matched_idx = moving_matched_idx(fixed_matched_idx);
        selected_dist = selected_dist(fixed_matched_idx);
    case {'mutual_closest'}
        [fixed_matched_idx, moving_matched_idx, selected_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);
end

registered_str.Fixed_sub = pc_fix_sub(fixed_matched_idx,:);
registered_str.Moving_sub_0 = pc_mov_sub(moving_matched_idx, :);
registered_str.Moving_sub = registered_str.Y(moving_matched_idx, :);
registered_str.Fixed_matched_idx = fixed_matched_idx;
registered_str.Moving_matched_idx = moving_matched_idx;
registered_str.Dist_fixed_2_moving = selected_dist; 
registered_str.Mutual_matched_rate = numel(fixed_matched_idx) ./ size(registered_str.Y, 1);
registered_str.Number_of_pair = numel(fixed_matched_idx);
registered_str = rmfield(registered_str, 'Y');
registered_str.establish_correspondance_method = opt.establish_correspondance_method;
end