function anisotropy_stat = fun_analysis_get_anisotropy_stat_from_vessel_graph(vessel_graph, ...
    selected_r_range, weight_method)

if nargin < 2
    selected_r_range = [0, inf];
    weight_method = 'none';
elseif nargin < 3
    weight_method = 'none';
end
volume_weighted_num_simulation = 10000;

if selected_r_range(1) > selected_r_range(2)
    error('The first element of selected_r_range should be smaller than the second one');
end

persistent uni_ori_isotropy_info str_template
if strcmp(weight_method, 'none') && isempty(uni_ori_isotropy_info)
    uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
end
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_data', 'cov_mat', ...
        'svd_u', 'svd_value', 'svd_min2max', 'svd_value_sum', 'svd_value_norm', ...
        'svd_value_ratio', 'svd_max_vec', 'min2max_z', 'svd_1_z', 'fa_z', ...
        'min2max_p', 'svd_1_p', 'fa_p', 'weight_method', 'corr_zr', 'corr_zr_p'});
end
anisotropy_stat = str_template;
anisotropy_stat.weight_method = weight_method;
% Compute the link features on the fly
num_l = vessel_graph.link.num_cc;
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
    if vessel_graph.link.num_voxel_per_cc(iter_link) > 1
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
end
lf.ep1_to_ep2_direction_vec = lf.ep1_to_ep2_direction_vec.';

selected_Q = lf.dt_median >= selected_r_range(1) & ...
    lf.dt_median <= selected_r_range(2);
switch weight_method
    case 'none'
        ep2ep_vec = lf.ep1_to_ep2_direction_vec(selected_Q, :);
    case 'volume'
        weight_vec = lf.length(selected_Q) .* (lf.dt_median(selected_Q) .^ 2); % Doesn't need to multiple pi here actually, since it won't effect the result. 
        ep2ep_vec = bsxfun(@times, lf.ep1_to_ep2_direction_vec(selected_Q, :), ...
            weight_vec);
    case 'length'
        weight_vec = lf.length(selected_Q);
        ep2ep_vec = bsxfun(@times, lf.ep1_to_ep2_direction_vec(selected_Q, :), ...
            weight_vec);
    otherwise
        error('Unknown weighted method');    
end

if nnz(selected_Q) < 3
    return;
end
anisotropy_stat = fun_analysis_get_link_anisotropy(ep2ep_vec, true);
anisotropy_stat.weight_method = weight_method;
% Check the correlation between selected link radius and the z-component of
% the orientation vector
[anisotropy_stat.corr_zr, anisotropy_stat.corr_zr_p] = corr(lf.dt_median(selected_Q), ...
    abs(lf.ep1_to_ep2_direction_vec(selected_Q, 3)));
if ~isempty(anisotropy_stat.svd_min2max)
    if strcmp(weight_method, 'none')
        anisotropy_stat.null.min2max_mean = uni_ori_isotropy_info.min2max.mean(anisotropy_stat.num_data);
        anisotropy_stat.null.min2max_std = uni_ori_isotropy_info.min2max.std(anisotropy_stat.num_data);
        
        anisotropy_stat.null.svd_1_mean = uni_ori_isotropy_info.svd_1.mean(anisotropy_stat.num_data);
        anisotropy_stat.null.svd_1_std = uni_ori_isotropy_info.svd_1.std(anisotropy_stat.num_data);
        
        anisotropy_stat.null.fa_mean = uni_ori_isotropy_info.fa.mean(anisotropy_stat.num_data);
        anisotropy_stat.null.fa_std = uni_ori_isotropy_info.fa.std(anisotropy_stat.num_data);
        
        tmp_simu_data.svd_min2max = [];
        tmp_simu_data.svd_1 = [];
        tmp_simu_data.fa = [];
    else
        tmp_uni_ori_str = fun_simulation_get_weighted_isotropy_stat(weight_vec, volume_weighted_num_simulation);
        anisotropy_stat.null.min2max_mean = tmp_uni_ori_str.svd_min2max_mean;
        anisotropy_stat.null.min2max_std = tmp_uni_ori_str.svd_min2max_std;
        
        anisotropy_stat.null.svd_1_mean = tmp_uni_ori_str.svd_ratio_mean(1);
        anisotropy_stat.null.svd_1_std = tmp_uni_ori_str.svd_ratio_std(1);
        
        anisotropy_stat.null.fa_mean = tmp_uni_ori_str.fractional_anisotropy_mean;
        anisotropy_stat.null.fa_std = tmp_uni_ori_str.fractional_anisotropy_std;
        
        tmp_simu_data = tmp_uni_ori_str.data;
    end
    anisotropy_stat.min2max_z = (anisotropy_stat.svd_min2max - anisotropy_stat.null.min2max_mean) / anisotropy_stat.null.min2max_std;
    anisotropy_stat.svd_1_z = (anisotropy_stat.svd_value_ratio(1) - anisotropy_stat.null.svd_1_mean) / anisotropy_stat.null.svd_1_std;
    anisotropy_stat.fa_z = (anisotropy_stat.fractional_anisotropy - anisotropy_stat.null.fa_mean) / anisotropy_stat.null.fa_std;
    % Compute p-value
    anisotropy_stat.min2max_p = nnz(anisotropy_stat.svd_min2max >= ...
        tmp_simu_data.svd_min2max) / numel(tmp_simu_data.svd_min2max);
    anisotropy_stat.svd_1_p = nnz(anisotropy_stat.svd_value_ratio(1) <= ...
        tmp_simu_data.svd_1) / numel(tmp_simu_data.svd_1);
    anisotropy_stat.fa_p = nnz(anisotropy_stat.fractional_anisotropy <= ...
        tmp_simu_data.fa) / numel(tmp_simu_data.fa);
end
%% Test the correlation between the z-component and the radius
% tmp_r = lf.dt_median;
% tmp_z = lf.ep1_to_ep2_direction_vec(:, 3);
% tmp_valid_Q = ~isnan(tmp_z) & ~isnan(tmp_r);
% tmp_corr = corr(tmp_r(tmp_valid_Q), tmp_z(tmp_valid_Q));
% fprintf('The correlation between the z-component and the radius is %f\n', tmp_corr);
end