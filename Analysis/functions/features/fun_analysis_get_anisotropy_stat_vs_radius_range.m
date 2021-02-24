function anisotropy_stat = fun_analysis_get_anisotropy_stat_vs_radius_range(vessel_graph, ...
    selected_r_range, weight_method)

volume_weighted_num_simulation = 10000;

persistent str_template
% persistent uni_ori_isotropy_info
% if strcmp(weight_method, 'none') && isempty(uni_ori_isotropy_info)
%     uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
% end
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_data', ...
        'svd_u', 'svd_value', 'svd_min2max', 'svd_value_sum', ...
        'svd_value_ratio', 'svd_max_vec', 'fa', 'svd_1_z', 'min2max_z', 'fa_z', ...
        'min2max_p', 'svd_1_p', 'fa_p',...
        'weight_method', 'weight_sum', 'select_r_min', 'select_r_max'});
end
anisotropy_stat = str_template;
anisotropy_stat.weight_method = weight_method;
% Compute the link features on the fly
is_not_single_voxel_link_Q = vessel_graph.link.num_voxel_per_cc > 1;
num_l = nnz(is_not_single_voxel_link_Q);
vessel_graph.link.cc_ind = vessel_graph.link.cc_ind(is_not_single_voxel_link_Q);
if num_l < 3
    return;
end

lf.ep1_to_ep2_direction_vec = nan(3, num_l);
lf.dt_median = nan(num_l, 1);

if ~strcmp(weight_method, 'none')
    compute_length_Q = true;
    lf.length = nan(num_l, 1);
else
    compute_length_Q = false;
end
for iter_link = 1 : num_l
        tmp_ind = vessel_graph.link.cc_ind{iter_link};
        tmp_sub = fun_ind2sub(vessel_graph.num.mask_size, tmp_ind);
        tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
        tmp_ep1_to_ep2_vec_norm = sqrt(sum((tmp_ep1_to_ep2_voxel_vec ).^2));
        lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec ./ ...
            tmp_ep1_to_ep2_vec_norm;
        tmp_dt = full(vessel_graph.radius(tmp_ind));
        lf.dt_median(iter_link) = median(tmp_dt(tmp_dt>0));
        if compute_length_Q
            lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);
        end
end
lf.ep1_to_ep2_direction_vec = lf.ep1_to_ep2_direction_vec.';
% Generate masking matrix, where each column correspond to each radius
% range
num_radius_range = size(selected_r_range, 2);
selected_mask = false(num_l, num_radius_range);
for iter_range = 1 : num_radius_range
    selected_mask(:, iter_range) = lf.dt_median >= selected_r_range(1, iter_range) ...
        & lf.dt_median <= selected_r_range(2, iter_range);
end
selected_overall = any(selected_mask, 2);
if nnz(selected_overall) <= 3
    return;
elseif ~all(selected_overall)
    lf = fun_structure_field_indexing(lf, selected_overall);
    selected_mask = selected_mask(selected_overall, :);
end

switch weight_method
    case 'none'
        weight_vec = ones(num_l, 1);
    case 'volume'
        weight_vec = lf.length .* (lf.dt_median .^ 2); % Doesn't need to multiple pi here actually, since it won't effect the result. 
    case 'length'
        weight_vec = lf.length;
    otherwise
        error('Unknown weighted method');    
end
ep2ep_vec = bsxfun(@times, lf.ep1_to_ep2_direction_vec, ...
            weight_vec);
weight_sum = sum(weight_vec .* selected_mask);

anisotropy_stat = fun_analysis_get_link_anisotropy_masking(ep2ep_vec, selected_mask);
anisotropy_stat.weight_sum = weight_sum;
anisotropy_stat.weight_method = weight_method;
anisotropy_stat.select_r_min = selected_r_range(1, :);
anisotropy_stat.select_r_max = selected_r_range(2, :);
if ~isempty(anisotropy_stat.svd_min2max)
    if strcmp(weight_method, 'none')
        error('Legacy option');
    else
        tmp_uni_ori_str = fun_simulation_get_weighted_isotropy_stat_with_masking(...
            weight_vec, selected_mask, volume_weighted_num_simulation);
        
        anisotropy_stat.null.svd_min2max = tmp_uni_ori_str.svd_min2max;
        
        anisotropy_stat.null.svd_1 = tmp_uni_ori_str.svd_1;
        
        anisotropy_stat.null.fa = tmp_uni_ori_str.fa;
                
        tmp_simu_data = tmp_uni_ori_str.data;
        tmp_simu_data.is_valid_Q = tmp_uni_ori_str.is_valid_Q;
    end
    anisotropy_stat.svd_min2max_z = (anisotropy_stat.svd_min2max - anisotropy_stat.null.svd_min2max.mean)...
        ./ anisotropy_stat.null.svd_min2max.std;
    anisotropy_stat.svd_1_z = (anisotropy_stat.svd_1 - anisotropy_stat.null.svd_1.mean)...
        ./ anisotropy_stat.null.svd_1.std;
    anisotropy_stat.fa_z = (anisotropy_stat.fa - anisotropy_stat.null.fa.mean)...
        ./ anisotropy_stat.null.fa.std;
    % Compute p-value
    anisotropy_stat.svd_min2max_p = sum(bsxfun(@ge, anisotropy_stat.svd_min2max, ...
        tmp_simu_data.svd_min2max), 1) ./ size(tmp_simu_data.svd_min2max, 1);
    anisotropy_stat.svd_1_p = sum(bsxfun(@le, anisotropy_stat.svd_1(1, :), ...
        tmp_simu_data.svd_1), 1) ./ size(tmp_simu_data.svd_1, 1);
    anisotropy_stat.fa_p = sum(bsxfun(@le, anisotropy_stat.fa, ...
        tmp_simu_data.fa), 1) ./ size(tmp_simu_data.fa, 1);
    
    if ~all(tmp_simu_data.is_valid_Q)
        anisotropy_stat.svd_min2max_p(~tmp_simu_data.is_valid_Q) = nan;
        anisotropy_stat.svd_1_p(~tmp_simu_data.is_valid_Q) = nan;
        anisotropy_stat.fa_p(~tmp_simu_data.is_valid_Q) = nan;
    end    
end
%% Test the correlation between the z-component and the radius
% tmp_r = lf.dt_median;
% tmp_z = lf.ep1_to_ep2_direction_vec(:, 3);
% tmp_valid_Q = ~isnan(tmp_z) & ~isnan(tmp_r);
% tmp_corr = corr(tmp_r(tmp_valid_Q), tmp_z(tmp_valid_Q));
% fprintf('The correlation between the z-component and the radius is %f\n', tmp_corr);
end