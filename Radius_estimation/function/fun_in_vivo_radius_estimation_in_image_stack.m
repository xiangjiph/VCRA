function est_radius = fun_in_vivo_radius_estimation_in_image_stack(invivo_im_prof_str, ...
    tile_image, image_voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um, est_method)

%% Parameters
max_r_to_est_um = max(invivo_im_prof_str.vsl_rds_list_um);
ori_num_voxel_half = 2;
ori_num_voxel = ori_num_voxel_half * 2 + 1;
%% Parse input
num_link = numel(skeleton_cc_sub_um);
assert(num_link == numel(skeleton_cc_r_um), 'Mismatch array size');
tile_im_size = size(tile_image);
%% Initialization
[link_ori_vec_z, best_fit_para, best_fit_corr, refine_r_est] = deal(cell(num_link, 1));
est_radius_stat = struct;
[est_radius_stat.median, est_radius_stat.mean, est_radius_stat.std] = deal(nan(4, num_link));
%% Computation
tmp_re_est_tic = tic;
for iter_link = 1 : num_link
    tmp_cc_sub_um = skeleton_cc_sub_um{iter_link};
    tmp_cc_r = double(skeleton_cc_r_um{iter_link});
    tmp_cc_num_voxel = numel(tmp_cc_r);
    assert(size(tmp_cc_sub_um, 1) == tmp_cc_num_voxel, 'Mismatch array size');
    tmp_cc_r_med = max(1, round(median(tmp_cc_r)));
%     tmp_cc_half_step = ceil(tmp_cc_r_med / 2 );
    tmp_cc_half_step = 0;
    if tmp_cc_r_med <= max_r_to_est_um
        % Determine the broadcasting indices
        if tmp_cc_num_voxel >= ori_num_voxel
            test_cc_r_est_idx = (ori_num_voxel_half + 1) : (2 * tmp_cc_half_step + 1) : (tmp_cc_num_voxel - ori_num_voxel_half);
            test_cc_r_seg_upper_ind = test_cc_r_est_idx + tmp_cc_half_step;
            test_cc_r_seg_upper_ind(end) = tmp_cc_num_voxel;
            test_cc_r_seg_num_vxl = diff([0, test_cc_r_seg_upper_ind]);
            assert(sum(test_cc_r_seg_num_vxl) == tmp_cc_num_voxel, 'The total length of the segments does not equal the length of the vessel');
        else
            test_cc_r_est_idx = ceil(tmp_cc_num_voxel / 2);
            test_cc_r_seg_num_vxl = tmp_cc_num_voxel;
        end
        num_r_est = numel(test_cc_r_est_idx);
        
        [tmp_cc_est_r_list] = deal(nan(3, num_r_est));
        [tmp_fit_corr, tmp_cc_ori_vec_z] = deal(nan(num_r_est, 1));
        [tmp_best_fit_parameter] = deal(nan(2, num_r_est));
        
        for iter_est = 1 : num_r_est
            tmp_voxel_idx = test_cc_r_est_idx(iter_est);
            tmp_voxel_r = tmp_cc_r(tmp_voxel_idx);
            % Coordinate of neighbor voxels - for estimating the tilt angle
            tmp_sample_ind_min = max(1, tmp_voxel_idx - ori_num_voxel_half);
            tmp_sample_ind_max = min(tmp_cc_num_voxel, tmp_voxel_idx + ori_num_voxel_half);
            tmp_neighbor_voxel_sub_um = tmp_cc_sub_um(tmp_sample_ind_min : tmp_sample_ind_max, :);
            % Compute the orientation vector
            tmp_vec = fun_radius_estimation_get_segment_orientation_vector(tmp_neighbor_voxel_sub_um);
            % Find the corresponding point in the full resolution image tile
            tmp_neighbor_voxel_sub_r = round(tmp_neighbor_voxel_sub_um ./ image_voxel_size_um);
            tmp_neighbor_voxel_sub_r = max(1, min(tile_im_size, tmp_neighbor_voxel_sub_r));
            
            para_est_str = fun_in_vivo_radius_parameters_joint_estimation(invivo_im_prof_str, tile_image, ...
                tmp_neighbor_voxel_sub_r, tmp_voxel_r, tmp_vec, ...
                image_voxel_size_um, est_method);
            
            if ~isempty(para_est_str)
                % Do not reject the estimation result. Results can be
                % post-selected later
%                 if abs(para_est_str.vsl_r_est - tmp_voxel_r) < max_est_update_diff_um
                tmp_cc_est_r_list(:, iter_est) = para_est_str.vsl_r_est;
                tmp_best_fit_parameter(:, iter_est) = para_est_str.best_fit_parameter;
                tmp_fit_corr(iter_est) = para_est_str.corr_max;
                tmp_cc_ori_vec_z(iter_est) = abs(tmp_vec(3));
%                 end
            end
        end
        %% Duplicate the radius estimation for connected component voxels
        tmp_cc_est_r_list = tmp_cc_est_r_list.';
        tmp_best_fit_parameter = tmp_best_fit_parameter.';

        % Use the estimated radius in selected skeleton voxels as the
        % radius of its neighboring skeleton voxels - duplicating
        tmp_cc_est_r_list = cat(2, tmp_cc_r, repelem(tmp_cc_est_r_list, test_cc_r_seg_num_vxl, 1));
        tmp_cc_ori_vec_z = repelem(tmp_cc_ori_vec_z, test_cc_r_seg_num_vxl, 1);
        tmp_fit_corr = repelem(tmp_fit_corr, test_cc_r_seg_num_vxl, 1);
        tmp_best_fit_parameter = repelem(tmp_best_fit_parameter, test_cc_r_seg_num_vxl, 1);
    else
        tmp_cc_est_r_list = cat(2, tmp_cc_r, nan(tmp_cc_num_voxel, 3));
        [tmp_fit_corr, tmp_cc_ori_vec_z] = deal(nan(tmp_cc_num_voxel, 1));
        tmp_best_fit_parameter = nan(tmp_cc_num_voxel, 2);
    end 
    %% Save cc estimation resoult
    link_ori_vec_z{iter_link} = tmp_cc_ori_vec_z;
    refine_r_est{iter_link} = tmp_cc_est_r_list;
    best_fit_para{iter_link} = tmp_best_fit_parameter;
    best_fit_corr{iter_link} = tmp_fit_corr;
    
%     fprintf('Finish estimating parameters near link %d\n', iter_link);
    % Evaluate the reliability of radius estimation before saving
    est_radius_stat.median(:, iter_link) = median(tmp_cc_est_r_list, 1, 'omitnan');
    est_radius_stat.mean(:, iter_link) = mean(tmp_cc_est_r_list, 1, 'omitnan'); 
    est_radius_stat.std(:, iter_link) = std(tmp_cc_est_r_list, 0, 1, 'omitnan');
end
fprintf('Finish estimation. Elapsed time is %f seconds.\n', toc(tmp_re_est_tic));
%% Return data
est_radius = struct;
est_radius.r = refine_r_est;
est_radius.local_ori_vec_z = link_ori_vec_z;
est_radius.best_fit_parameters = best_fit_para;
est_radius.best_fit_corr = best_fit_corr;

est_radius_stat.median = est_radius_stat.median.';
est_radius_stat.mean = est_radius_stat.mean.';
est_radius_stat.std = est_radius_stat.std.';

est_radius.cc_stat = est_radius_stat;
end