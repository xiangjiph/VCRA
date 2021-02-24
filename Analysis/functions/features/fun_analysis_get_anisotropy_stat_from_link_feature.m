function anisotropy_stat = fun_analysis_get_anisotropy_stat_from_link_feature(link_feature, max_select_radius, vol_weighted_Q)

if nargin < 2
    max_select_radius = inf;
    vol_weighted_Q = false;
elseif nargin < 3
    vol_weighted_Q = false;
end
volume_weighted_num_simulation = 10000; % For normal distribution, the average and std has uncertainty at about 1%.

persistent uni_ori_isotropy_info str_template
if ~vol_weighted_Q && isempty(uni_ori_isotropy_info)
    uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
end

if ~isempty(link_feature)
    if isfinite(max_select_radius)
        tmp_radius = link_feature.dt_median;
        tmp_selected_Q = tmp_radius <= max_select_radius;
    else
        tmp_selected_Q = true(size(link_feature.dt_median));
    end
    
    if vol_weighted_Q
        tmp_ep2ep_vec = bsxfun(@times, link_feature.ep1_to_ep2_direction_vec, link_feature.in_bbox_volume);
    else
        tmp_ep2ep_vec = link_feature.ep1_to_ep2_direction_vec;
    end
    
    anisotropy_stat = fun_analysis_get_link_anisotropy(tmp_ep2ep_vec(tmp_selected_Q, :), true);
    if ~isempty(anisotropy_stat.svd_min2max)
        if ~vol_weighted_Q
            anisotropy_stat.null.min2max_mean = uni_ori_isotropy_info.min2max.mean(anisotropy_stat.num_data);
            anisotropy_stat.null.min2max_std = uni_ori_isotropy_info.min2max.std(anisotropy_stat.num_data);
            
            anisotropy_stat.null.svd_1_mean = uni_ori_isotropy_info.svd_1.mean(anisotropy_stat.num_data);
            anisotropy_stat.null.svd_1_std = uni_ori_isotropy_info.svd_1.std(anisotropy_stat.num_data);
            
            anisotropy_stat.null.fa_mean = uni_ori_isotropy_info.fa.mean(anisotropy_stat.num_data);
            anisotropy_stat.null.fa_std = uni_ori_isotropy_info.fa.std(anisotropy_stat.num_data);
        else
            tmp_uni_ori_str = fun_simulation_get_weighted_isotropy_stat(lf.volume, volume_weighted_num_simulation);
            anisotropy_stat.null.min2max_mean = tmp_uni_ori_str.svd_min2max_mean;
            anisotropy_stat.null.min2max_std = tmp_uni_ori_str.svd_min2max_std;
            
            anisotropy_stat.null.svd_1_mean = tmp_uni_ori_str.svd_ratio_mean(1);
            anisotropy_stat.null.svd_1_std = tmp_uni_ori_str.svd_ratio_std(1);
            
            anisotropy_stat.null.fa_mean = tmp_uni_ori_str.fractional_anisotropy_mean;
            anisotropy_stat.null.fa_std = tmp_uni_ori_str.fractional_anisotropy_std;
        end
        anisotropy_stat.min2max_z = (anisotropy_stat.svd_min2max - anisotropy_stat.null.min2max_mean) / anisotropy_stat.null.min2max_std;
        anisotropy_stat.svd_1_z = (anisotropy_stat.svd_value_ratio(1) - anisotropy_stat.null.svd_1_mean) / anisotropy_stat.null.svd_1_std;
        anisotropy_stat.fa_z = (anisotropy_stat.fractional_anisotropy - anisotropy_stat.null.fa_mean) / anisotropy_stat.null.fa_std;
    end
else
    if isempty(str_template)
        str_template = fun_initialized_structure_array_with_fieldname_list({'num_data', 'cov_mat', ...
            'svd_u', 'svd_value', 'svd_min2max', 'svd_value_sum', 'svd_value_norm', ...
            'svd_value_ratio', 'svd_max_vec', 'min2max_z', 'svd_1_z'});
    end
    anisotropy_stat = str_template;    
end
end