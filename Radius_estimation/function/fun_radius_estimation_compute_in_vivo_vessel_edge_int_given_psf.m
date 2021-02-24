function edge_int_str = fun_radius_estimation_compute_in_vivo_vessel_edge_int_given_psf(psf_sigma_um, vsl_r_list_um, vsl_theta_list, b_lyr_thickness_um, rbc_int_n_list, dx_um, need_vertical_prof_Q)

if nargin < 7
    need_vertical_prof_Q = true;
end
num_rbc_int_n = numel(rbc_int_n_list);
num_theta = numel(vsl_theta_list);
num_vsl_r = numel(vsl_r_list_um);
psf_FWHM_um = psf_sigma_um .* (2*sqrt(2*log(2)));
dx_list = nan(size(vsl_r_list_um));
[abs_min_edge_int, abs_center_int, abs_max_int, center_normalized_min_edge_int, max_normalized_min_edge_int] = deal(nan(num_rbc_int_n, num_theta, num_vsl_r));
[max_n_radial_int_dist, max_n_vertical_radial_int_dist, ...
    center_n_vertical_radial_int_dist, center_n_radial_int_dist] = deal(cell(num_rbc_int_n, num_theta, num_vsl_r));
%% Calculation 
% Start from large value to check RAM overflow on GPU. 
for tmp_r_idx = num_vsl_r : -1 : 1 
    vessel_radius = vsl_r_list_um(tmp_r_idx);
    fprintf('Simulating image of vessel (r = %.2f um) with PSF FWHM (%.2f, %.2f, %.2f) um\n', vessel_radius, psf_FWHM_um);
    dx_list(tmp_r_idx) = dx_um;   
    vessel_radius_in_voxel = round(vessel_radius / dx_um);
    vessel_bc_t_in_voxel = round(b_lyr_thickness_um / dx_um);
    if vessel_bc_t_in_voxel >= vessel_radius_in_voxel || b_lyr_thickness_um >= vessel_radius
        continue;
    end
    
    psf_sigma_in_voxel = psf_sigma_um ./ dx_um;
    % Generate gaussian profile PSF
    psf = fun_gaussian_kernel(psf_sigma_in_voxel);
    ker1 = single(fun_gaussian_kernel(psf_sigma_in_voxel(1)));
    ker2 = single(fun_gaussian_kernel(psf_sigma_in_voxel(2))');
    ker3 = single(reshape(fun_gaussian_kernel(psf_sigma_in_voxel(3)),1,1,[]));
    
    psf_size = size(psf);
    tmp_tic = tic;
    for tmp_th_idx = 1 : num_theta
%         profile on 
        theta = vsl_theta_list(tmp_th_idx);
        if need_vertical_prof_Q
            % The entire PSF size is required in each half-size because the
            % object above the cutoff boundary also contributes to the
            % convolution value. Half-size can be used for x-y dimension,
            % but since the PSF is much smaller in xy direction, we also
            % use 
            simu_system_half_size = nan(1,3, 'single');
            simu_system_half_size(1:2) = ceil(vessel_radius_in_voxel + psf_size(1:2) + 1);
            z_extend_length_pxl = ceil( 4 / dx_um); % Get the extra profile for PSF fitting
            % the factor of 2 is chosen to get the valid intensity profile for 60\deg elevation angle
            simu_system_half_size(3) = ceil(2 * vessel_radius_in_voxel + z_extend_length_pxl + psf_size(3) + 1); 
        else
            simu_system_half_size = ones(1,3) * vessel_radius_in_voxel;
            simu_system_half_size(1) = simu_system_half_size(1) + psf_size(1);
            simu_system_half_size(2) = simu_system_half_size(2) + psf_size(2);
            simu_system_half_size(3) = psf_size(3) + 1;            
            simu_system_half_size = single(ceil(simu_system_half_size));
        end        
        vessel_im_prof_str = fun_radius_estimation_get_theoretical_in_vivo_vessel_im_profile(...
            vessel_radius_in_voxel, vessel_bc_t_in_voxel, ...
            simu_system_half_size, theta, ker1, ker2, ker3, true);
        vessel_center_position = vessel_im_prof_str.vessel_center_sub;
        for iter_rbc_n_int = 1 : num_rbc_int_n
        % Voxel subscripts on the xy plane that pass the center:
        % Extract the edge intensity from the simulaton image
        % Intensity distribution along the radius on orizontal plane
            tmp_int_subtraction_factor = 1 - rbc_int_n_list(iter_rbc_n_int);            
            tmp_vsl_w_rbc_image = vessel_im_prof_str.vessel_image - vessel_im_prof_str.vessel_rbc_tube_image * tmp_int_subtraction_factor;
            
            radial_intensity_distribution = tmp_vsl_w_rbc_image(vessel_center_position(1):end , vessel_center_position(2), vessel_center_position(3));
            % Intensity distribution along the radius on the vertical plane
            vertical_radial_intensity_distribution = squeeze(tmp_vsl_w_rbc_image(...
                vessel_center_position(1), vessel_center_position(2), ...
                vessel_center_position(3) : end));

            radial_pos = (0:(numel(radial_intensity_distribution) - 1))' .* dx_um;
            tmp_max_int = max(radial_intensity_distribution(:));
            abs_max_int(iter_rbc_n_int, tmp_th_idx, tmp_r_idx) = tmp_max_int;
            max_n_radial_int_dist{iter_rbc_n_int, tmp_th_idx, tmp_r_idx} = cat(2, radial_pos, radial_intensity_distribution ./ tmp_max_int);
            
            tmp_abs_center_int = tmp_vsl_w_rbc_image(vessel_center_position(1), vessel_center_position(2), vessel_center_position(3));
            assert(tmp_abs_center_int ==  radial_intensity_distribution(1));
            center_n_radial_int_dist{iter_rbc_n_int, tmp_th_idx, tmp_r_idx} = cat(2, radial_pos, radial_intensity_distribution ./ tmp_abs_center_int);

            vertical_radial_pos = (0 : numel(vertical_radial_intensity_distribution) - 1)' .* dx_um;
            center_n_vertical_radial_int_dist{iter_rbc_n_int, tmp_th_idx, tmp_r_idx} = cat(2, vertical_radial_pos, ...
                vertical_radial_intensity_distribution ./ vertical_radial_intensity_distribution(1));
            max_n_vertical_radial_int_dist{iter_rbc_n_int, tmp_th_idx, tmp_r_idx} = cat(2, vertical_radial_pos, ...
                vertical_radial_intensity_distribution ./ max(vertical_radial_intensity_distribution(:)));

            assert(abs(radial_pos(vessel_radius_in_voxel + 1) - vessel_radius) < 1e-3); % The simulation should generate vessel of size excatly equals to vessel_raidus
            tmp_abs_lateral_edge_int = radial_intensity_distribution(vessel_radius_in_voxel+1);
            abs_min_edge_int(iter_rbc_n_int, tmp_th_idx, tmp_r_idx) = tmp_abs_lateral_edge_int;
            
            abs_center_int(iter_rbc_n_int, tmp_th_idx, tmp_r_idx) = tmp_abs_center_int;
            center_normalized_min_edge_int(iter_rbc_n_int, tmp_th_idx, tmp_r_idx) = tmp_abs_lateral_edge_int / tmp_abs_center_int;
            max_normalized_min_edge_int(iter_rbc_n_int, tmp_th_idx, tmp_r_idx) = tmp_abs_lateral_edge_int / tmp_max_int;
        end
%         profile off
%         profile viewer;
        %% Visualization
%         tmp_plot_data = cat(2, center_n_radial_int_dist{:, tmp_th_idx, tmp_r_idx});
%         tmp_y = tmp_plot_data(:, 2:2:end);
%         plot(tmp_y);
%         legend;
    end
    fprintf('Elapsed time is %f seconds\n', toc(tmp_tic));
end
%% Save result
edge_int_str = struct;
edge_int_str.psf_sigma_um = psf_sigma_um;
edge_int_str.psf_FWHM_um = psf_FWHM_um;
edge_int_str.vsl_boundary_layer_thickness_um = b_lyr_thickness_um;
edge_int_str.vsl_rbc_int_n = rbc_int_n_list;
edge_int_str.vsl_r_list = vsl_r_list_um;
edge_int_str.vsl_theta_list = vsl_theta_list;
edge_int_str.psf_sigma_um = psf_sigma_um;
edge_int_str.simulation_dx_min_scale_ratio = dx_um;
edge_int_str.simulation_voxel_size_list = dx_list;
edge_int_str.abs_min_edge_int = abs_min_edge_int;
edge_int_str.abs_center_int = abs_center_int;
edge_int_str.normalized_min_edge_int = center_normalized_min_edge_int;
edge_int_str.normalized_radial_int_dist = center_n_radial_int_dist;
edge_int_str.normalized_vertical_radial_int_dist = center_n_vertical_radial_int_dist;

edge_int_str.max_normalized_min_edge_int = max_normalized_min_edge_int;
edge_int_str.max_normalized_radial_int_dist = max_n_radial_int_dist;
edge_int_str.max_normalized_vertical_radial_int_dist = max_n_vertical_radial_int_dist;

% tmp_int = struct;
% tmp_z_comp = sin(vsl_theta_list);
% [tmp_int.z_comp, tmp_int.radius] = ndgrid(tmp_z_comp, vsl_r_list_um);
% tmp_int.variable_name = {'z_component', 'radius'};
% tmp_int.Interpolation_method = 'linear';
% tmp_int.Extrapolation_method = 'nearest';
% tmp_int.n_min_edge_int = griddedInterpolant(tmp_int.z_comp, tmp_int.radius, edge_int_str.normalized_min_edge_int, ...
%     tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
% edge_int_str.n_min_edge_int_interpolation = tmp_int;
% 
% tmp_int.n_min_edge_int = griddedInterpolant(tmp_int.z_comp, tmp_int.radius, edge_int_str.max_normalized_min_edge_int, ...
%     tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
% edge_int_str.max_n_min_edge_int_interpolation = tmp_int;
end