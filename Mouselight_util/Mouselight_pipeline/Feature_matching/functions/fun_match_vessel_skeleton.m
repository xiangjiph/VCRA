function registered_str = fun_match_vessel_skeleton(...
    skel_fix_sub, skel_moving_sub, skel_fix_label, skel_moving_label, disp_vec_0, iadj, matchparams)
%% Documentation


%% Default setting
if isfield(matchparams, 'max_num_desc')
    sample_opt.total_num_descriptor = matchparams.max_num_desc;
else
    sample_opt.total_num_descriptor = 5000;
end
if ~isfield(matchparams, 'viz')
    vis_Q = false;
else
    vis_Q = matchparams.viz;
end
if isfield(matchparams, 'remove_descriptor_by_pdist2')
    sample_opt.rm_descriptor_by_pdist2_Q = matchparams.remove_descriptor_by_pdist2;
else
    sample_opt.rm_descriptor_by_pdist2_Q = true;
end
if ~isfield(matchparams, 'pixshift_search_shift_z')
    pixshift_search_shift_z = [1; -1] * [15, 30, 50];
else
    pixshift_search_shift_z = [1; -1] * matchparams.pixshift_search_shift_z;
end
switch iadj
    case 1
        sample_opt.max_disp_pixel_yxz = [15, 10, 5];
    case 2
        sample_opt.max_disp_pixel_yxz = [10, 15, 5];
    case 3
        sample_opt.max_disp_pixel_yxz = [30, 30, 30]; % I am not sure if these numbers are good 
end

pixshift_search_record = cell(size(pixshift_search_shift_z));
pixshift_search_num_pair = zeros(size(pixshift_search_shift_z));
num_pixel_search = numel(pixshift_search_record);
% Inconsistancy threshold
th_inconsistancy = 0.1;
% Parse input
registered_str = [];
disp_vec_0 = round(disp_vec_0);
disp_vec = disp_vec_0;

num_search_option = numel(pixshift_search_shift_z);
flag_stop = false;
iter = 0;

R_consistant = zeros(1,50);
exhaust_pixel_search_idx = 0;
while ~flag_stop && iter <= num_search_option % run a search
    iter = iter + 1;
%% Shift the pixels according to the stage displacement 
    des_1_sub = skel_fix_sub;
    des_1_label = skel_fix_label;
    
    des_1_sub_min = min(des_1_sub, [], 1);
    des_1_sub_max = max(des_1_sub, [], 1);
    
    des_2_sub_shift = bsxfun(@plus, skel_moving_sub, round(disp_vec));
    des_2_label = skel_moving_label;
    
    des_2_sub_shift_min = min(des_2_sub_shift, [], 1);
    des_2_sub_shift_max = max(des_2_sub_shift, [], 1);
    
    overlap_bbox_max = min([des_1_sub_max; des_2_sub_shift_max], [], 1);
    overlap_bbox_min = max([des_1_sub_min; des_2_sub_shift_min; [1,1,1]],[], 1);
    
    desc_1_selected_idx = find(all(bsxfun(@ge, des_1_sub, overlap_bbox_min) & bsxfun(@le, des_1_sub, overlap_bbox_max), 2));
    desc_2_selected_idx = find(all(bsxfun(@ge, des_2_sub_shift, overlap_bbox_min) & bsxfun(@le, des_2_sub_shift, overlap_bbox_max), 2));
    
    des_1_sub = des_1_sub(desc_1_selected_idx, :);
    des_2_sub_shift = des_2_sub_shift(desc_2_selected_idx, :);    
    if isempty(des_1_sub) || isempty(des_2_sub_shift)
        return;
    end
%% Sample the voxels for matching
    [des_1_sub, des_2_sub_shift, des_1_list_ind, des_2_list_ind] = fun_feature_match_sample_skeleton_voxel(des_1_sub, ...
        des_2_sub_shift, des_1_label, des_2_label, sample_opt);
    desc_1_selected_idx = desc_1_selected_idx(des_1_list_ind);
    desc_2_selected_idx = desc_2_selected_idx(des_2_list_ind);
    
    des_1_label = des_1_label(desc_1_selected_idx);
    des_2_label = des_2_label(desc_2_selected_idx);            
%% Coherent Point drift
    if size(des_1_sub,1)<3 || size(des_2_sub_shift,1)<3% not enough sample to match
        flag_stop = 1;
    else
%% Match the descriptor - nonregid        
        registered_str = fun_match_point_cloud_by_cpd(des_1_sub, ...
            des_2_sub_shift, matchparams);
        delete(gcf);
        mutual_matched_rate = registered_str.Mutual_matched_rate;
        num_matched_point = numel(registered_str.Fixed_matched_idx);
        if num_matched_point == 0
            return;
        end
        %% Outlier rejections - accoding to connected component
        % Check if the matched points are from the same connected
        % components or not. A good matches should contains many voxel
        % pairs from single connected components        
        des_1_label = des_1_label(registered_str.Fixed_matched_idx);
        des_2_label = des_2_label(registered_str.Moving_matched_idx);
        % Connected components in the input point set
        X_cc_idx = fun_bin_data_to_idx_list(des_1_label);
        X_cc_num = numel(X_cc_idx);
        X_cc_size = cellfun(@numel, X_cc_idx);
        % Define inconsistancy as the number of labels / number of voxels
        X_cc_inconsistancy = zeros(X_cc_num, 1);
        for iter1 = 1 : X_cc_num
            X_cc_inconsistancy(iter1) = numel(unique(des_2_label(X_cc_idx{iter1})))/X_cc_size(iter1);
        end
        matched_cc_x = X_cc_inconsistancy <= th_inconsistancy;

        Y_cc_idx = fun_bin_data_to_idx_list(des_2_label);
        Y_cc_num = numel(Y_cc_idx);
        Y_cc_size = cellfun(@numel, Y_cc_idx);
        % Define inconsistancy as the number of labels / number of voxels
        Y_cc_inconsistancy = zeros(Y_cc_num, 1);
        for iter1 = 1 : Y_cc_num
            Y_cc_inconsistancy(iter1) = numel(unique(des_2_label(Y_cc_idx{iter1})))/Y_cc_size(iter1);
        end
        matched_cc_y = Y_cc_inconsistancy <= th_inconsistancy;
        
        match_x_Q = false(num_matched_point, 1);
        match_x_Q(cat(2, X_cc_idx{matched_cc_x})) = true;
        
        match_y_Q = false(num_matched_point, 1);
        match_y_Q(cat(2, Y_cc_idx{matched_cc_y})) = true;
        
        matched_Q = match_x_Q & match_y_Q;
        %%
        tmp_num_matched_voxel = nnz(matched_Q);
        registered_str.Fixed_cc_label = des_1_label(matched_Q);
        registered_str.Moving_cc_label = des_2_label(matched_Q);
        
        registered_str = fun_structure_field_indexing(registered_str, matched_Q, ...
            {'Fixed_sub', 'Moving_sub_0', 'Moving_sub', 'Fixed_matched_idx', ...
            'Moving_matched_idx', 'Dist_fixed_2_moving'});
        
        if tmp_num_matched_voxel < 3
            consistent_rate = 0; % Too less points matched
        else
            consistent_rate = tmp_num_matched_voxel / numel(matched_Q);
        end
        R_consistant(iter) = consistent_rate;
        %% 
        % The displacement between the matched voxels should be consistant.
        % Remove the outlier pairs whose displacements are significant
        % different from the median displacement of the matched pairs
        disp_X_Y = registered_str.Fixed_sub - registered_str.Moving_sub_0 + disp_vec;
        disp_X_Y_med = median(disp_X_Y, 1);
        disp_X_Y_dev = bsxfun(@minus, disp_X_Y , disp_X_Y_med);
        disp_X_Y_dev_std = std(single(disp_X_Y_dev), 1);
        disp_X_Y_tol = min(15, max(5,disp_X_Y_dev_std * 3));
        disp_kept_Q = all(bsxfun(@le, abs(disp_X_Y_dev), disp_X_Y_tol), 2);
        
        registered_str = fun_structure_field_indexing(registered_str, disp_kept_Q, ...
            {'Fixed_sub', 'Moving_sub_0', 'Moving_sub', 'Fixed_matched_idx', ...
            'Moving_matched_idx', 'Dist_fixed_2_moving', ...
            'Fixed_cc_label', 'Moving_cc_label'});
        registered_str.Number_of_pair = nnz(disp_kept_Q);
        registered_str.Final_disp_vec = disp_X_Y_med;
        if ~isempty(registered_str.Fixed_sub) && size(registered_str.Fixed_sub, 1) > 100 %&& search_by_brute_force_Q
            search_by_brute_force_Q = false;
        else 
            search_by_brute_force_Q = true;
        end
%% Re-estiamte the initial pixel shift - Brute force       
        if matchparams.scan_pixshift_Q && mutual_matched_rate < 0.95
            disp('Matching not good enough. Shift the overlapping region and search for pairs again');
            if search_by_brute_force_Q
                exhaust_pixel_search_idx = exhaust_pixel_search_idx + 1;
                if exhaust_pixel_search_idx <= num_pixel_search
                    pixshift_search_record{exhaust_pixel_search_idx} = registered_str;
                    disp_vec = disp_vec_0;
                    disp_vec(iadj) = disp_vec(iadj) + pixshift_search_shift_z(exhaust_pixel_search_idx);
                    pixshift_search_num_pair(exhaust_pixel_search_idx) = nnz(disp_kept_Q);
                else
                    flag_stop = true;
                    [~, return_idx] = max(pixshift_search_num_pair);
                    registered_str = pixshift_search_record{return_idx};
                end
            else
                if iter == 1
                    registered_str_0 = registered_str;
                    disp_vec = disp_X_Y_med;
                else                  
                    if R_consistant(iter) > R_consistant(iter-1) && any(abs(disp_X_Y_med - disp_vec) > 5)
                        % If the new matching is better and is
                        % significantly different from the original initial
                        % displacement estimation.
                        registered_str_0 = registered_str;
                        disp_vec = disp_X_Y_med;
                    else
                        flag_stop = true;
                        registered_str = registered_str_0;
                    end
                end
            end
        else
            flag_stop = true;
            disp_vec = round(disp_X_Y_med);
        end
    end
end
registered_str.R_consistant = R_consistant;
end
