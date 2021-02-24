function edge_int_str = fun_radius_estimation_compute_vessel_edge_int_given_psf(psf_sigma_um, vsl_r_list, vsl_theta_list, simu_dx_min_scale_ratio, need_vertical_prof_Q)

if nargin < 5
    need_vertical_prof_Q = true;
end

psf_FWHM_um = psf_sigma_um .* (2*sqrt(2*log(2)));
simu_target_min_scale_div = round(1 ./ simu_dx_min_scale_ratio);
dx_list = nan(size(vsl_r_list));
[abs_min_edge_int, abs_center_int, normalized_min_edge_int] = deal(zeros(numel(vsl_theta_list), numel(vsl_r_list)));
[n_radial_int_dist, n_vertical_radial_int_dist] = deal(cell(numel(vsl_theta_list), numel(vsl_r_list)));
%% Calculation 
% Start from large value to check RAM overflow on GPU. 
for tmp_r_idx = numel(vsl_r_list) : -1 : 1 
    vessel_radius = vsl_r_list(tmp_r_idx);
    fprintf('Simulating image of vessel (r = %.2f um) with PSF FWHM (%.2f, %.2f, %.2f) um\n', vessel_radius, psf_FWHM_um);
    % Search for the simulation voxel size
    tmp_num_div = simu_target_min_scale_div;
    while true
        tmp_dx = vessel_radius / tmp_num_div;
        if all(psf_sigma_um .* 2 ./ tmp_dx >= simu_target_min_scale_div)
            fprintf('Set the voxle size for simulation to be %f um\n', tmp_dx);
            break;
        else
            tmp_num_div = tmp_num_div + 1;
        end
    end
    assert(abs(mod(vessel_radius, tmp_dx)) < 1e-4, 'Vessel radius / voxel size is not an integer');
    dx_list(tmp_r_idx) = tmp_dx;   
    
    psf_sigma_in_voxel = psf_sigma_um ./ tmp_dx;
    % Generate gaussian profile PSF
    psf = fun_gaussian_kernel(psf_sigma_in_voxel);
    ker1 = single(fun_gaussian_kernel(psf_sigma_in_voxel(1)));
    ker2 = single(fun_gaussian_kernel(psf_sigma_in_voxel(2))');
    ker3 = single(reshape(fun_gaussian_kernel(psf_sigma_in_voxel(3)),1,1,[]));
    
    psf_size = size(psf);
    tmp_tic = tic;
    for tmp_th_idx = 1 : numel(vsl_theta_list)
        theta = vsl_theta_list(tmp_th_idx);
        vessel_radius_in_voxel = round(vessel_radius / tmp_dx);        
        if need_vertical_prof_Q
            % The entire PSF size is required in each half-size because the
            % object above the cutoff boundary also contributes to the
            % convolution value. Half-size can be used for x-y dimension,
            % but since the PSF is much smaller in xy direction, we also
            % use 
            simu_system_half_size = nan(1,3, 'single');
            simu_system_half_size(1:2) = ceil(vessel_radius_in_voxel + psf_size(1:2) + 1);
            z_extend_length_pxl = ceil( 4 / tmp_dx); % Get the extra profile for PSF fitting
            % the factor of 2 is chosen to get the valid intensity profile for 60\deg elevation angle
            simu_system_half_size(3) = ceil(2 * vessel_radius_in_voxel + z_extend_length_pxl + psf_size(3) + 1); 
        else
            simu_system_half_size = ones(1,3) * vessel_radius_in_voxel;
            simu_system_half_size(1) = simu_system_half_size(1) + psf_size(1);
            simu_system_half_size(2) = simu_system_half_size(2) + psf_size(2);
            simu_system_half_size(3) = psf_size(3) + 1;            
            simu_system_half_size = single(ceil(simu_system_half_size));
        end        
        vessel_im_prof_str = fun_radius_estimation_get_theoretical_vessel_im_profile(...
            vessel_radius_in_voxel, simu_system_half_size, theta, ker1, ker2, ker3, true);
        
        % Voxel subscripts on the xy plane that pass the center:
        % Extract the edge intensity from the simulaton image
        % Intensity distribution along the radius on orizontal plane
        vessel_image = vessel_im_prof_str.vessel_image;        
        vessel_center_position = vessel_im_prof_str.vessel_center_sub;
        radial_intensity_distribution = vessel_image(vessel_center_position(1):end , vessel_center_position(2), vessel_center_position(3));
        % Intensity distribution along the radius on the vertical plane
        vertical_radial_intensity_distribution = squeeze(vessel_image(...
            vessel_center_position(1), vessel_center_position(2), ...
            vessel_center_position(3) : end));
        
        radial_pos = (0:(numel(radial_intensity_distribution) - 1))' .* tmp_dx;
        n_radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, radial_pos, radial_intensity_distribution ./ radial_intensity_distribution(1));
        
        vertical_radial_pos = (0 : numel(vertical_radial_intensity_distribution) - 1)' .* tmp_dx;
        n_vertical_radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, vertical_radial_pos, ...
            vertical_radial_intensity_distribution ./ vertical_radial_intensity_distribution(1));
        
        assert(abs(radial_pos(vessel_radius_in_voxel + 1) - vessel_radius) < 1e-3); % The simulation should generate vessel of size excatly equals to vessel_raidus 
        abs_lateral_edge_int = radial_intensity_distribution(vessel_radius_in_voxel+1);
        abs_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_lateral_edge_int;
        abs_center_int(tmp_th_idx, tmp_r_idx) = vessel_image(vessel_center_position(1), vessel_center_position(2), vessel_center_position(3));
        normalized_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_min_edge_int(tmp_th_idx,tmp_r_idx)/abs_center_int(tmp_th_idx, tmp_r_idx);
    end
    fprintf('Elapsed time is %f seconds\n', toc(tmp_tic));
end
%% Save result
edge_int_str = struct;
edge_int_str.psf_sigma_um = psf_sigma_um;
edge_int_str.psf_FWHM_um = psf_FWHM_um;
edge_int_str.vsl_r_list = vsl_r_list;
edge_int_str.vsl_theta_list = vsl_theta_list;
edge_int_str.psf_sigma_um = psf_sigma_um;
edge_int_str.simulation_dx_min_scale_ratio = simu_dx_min_scale_ratio;
edge_int_str.simulation_voxel_size_list = dx_list;
edge_int_str.abs_min_edge_int = abs_min_edge_int;
edge_int_str.abs_center_int = abs_center_int;
edge_int_str.normalized_min_edge_int = normalized_min_edge_int;
edge_int_str.normalized_radial_int_dist = n_radial_int_dist;
edge_int_str.normalized_vertical_radial_int_dist = n_vertical_radial_int_dist;

tmp_int = struct;
tmp_r = vsl_r_list;
tmp_z_comp = sin(vsl_theta_list);
[tmp_int.z_comp, tmp_int.radius] = ndgrid(tmp_z_comp, tmp_r);
tmp_int.variable_name = {'z_component', 'radius'};
tmp_int.Interpolation_method = 'linear';
tmp_int.Extrapolation_method = 'nearest';
tmp_int.n_min_edge_int = griddedInterpolant(tmp_int.z_comp, tmp_int.radius, edge_int_str.normalized_min_edge_int, ...
    tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
edge_int_str.n_min_edge_int_interpolation = tmp_int;
end