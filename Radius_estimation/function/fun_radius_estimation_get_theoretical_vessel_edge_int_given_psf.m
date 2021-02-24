function edge_int_str = fun_radius_estimation_get_theoretical_vessel_edge_int_given_psf(psf_FWHM_um, vsl_r_list, vsl_theta_list, dx)

if nargin < 4
    dx = min(psf_FWHM_um) / 10;
end

psf_sigma = psf_FWHM_um./(2*sqrt(2*log(2)));
psf_sigma_in_voxel = round(psf_sigma ./ dx);
% Generate gaussian profile PSF
psf = fun_gaussian_kernel(psf_sigma_in_voxel);
ker1 = fun_gaussian_kernel(psf_sigma_in_voxel(1));
ker2 = fun_gaussian_kernel(psf_sigma_in_voxel(2))';
ker3 = reshape(fun_gaussian_kernel(psf_sigma_in_voxel(3)),1,1,[]);

psf_size = size(psf);
[abs_min_edge_int, abs_center_int, normalized_min_edge_int] = deal(zeros(numel(vsl_theta_list), numel(vsl_r_list)));
[n_radial_int_dist, n_vertical_radial_int_dist] = deal(cell(numel(vsl_theta_list), numel(vsl_r_list)));
%% Calculation 
% Start from large value to check RAM overflow on GPU. 
for tmp_r_idx = numel(vsl_r_list) : -1 : 1 
    vessel_radius = vsl_r_list(tmp_r_idx);
    tmp_tic = tic;
    for tmp_th_idx = 1 : numel(vsl_theta_list)
        theta = vsl_theta_list(tmp_th_idx);
        vessel_radius_in_voxel = round(vessel_radius / dx);
        vessel_length_in_voxel = psf_size(3);        
        vessel_im_prof_str = fun_radius_estimation_get_theoretical_vessel_im_profile(...
            vessel_radius_in_voxel, vessel_length_in_voxel, theta, ker1, ker2, ker3, true);
        % Voxel subscripts on the xy plane that pass the center:
        % Extract the edge intensity from the simulaton image
        % Intensity distribution along the radius on orizontal plane
        vessel_image = vessel_im_prof_str.vessel_image;        
        vessel_center_position = vessel_im_prof_str.vessel_center_sub;
        lateral_radial_intensity_distribution = vessel_image(vessel_center_position(1):end , vessel_center_position(2), vessel_center_position(3));
        % Intensity distribution along the radius on the vertical plane
        vertical_radial_intensity_distribution = squeeze(vessel_image(...
            vessel_center_position(1), vessel_center_position(2), ...
            vessel_center_position(3) : end));
        
        lateral_radial_pos = (0:(numel(lateral_radial_intensity_distribution) - 1))' .* dx;
        lateral_int_profile_n = lateral_radial_intensity_distribution ./ lateral_radial_intensity_distribution(1);
        
        n_radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, lateral_radial_pos, lateral_int_profile_n);
        
        vertical_radial_pos = (0 : numel(vertical_radial_intensity_distribution) - 1)' .* dx;
        n_vertical_radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, vertical_radial_pos, ...
            vertical_radial_intensity_distribution ./ vertical_radial_intensity_distribution(1));
        
        abs_lateral_edge_int = lateral_radial_intensity_distribution(vessel_radius_in_voxel+1);
        abs_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_lateral_edge_int;
        abs_center_int(tmp_th_idx,tmp_r_idx) = vessel_image(vessel_center_position(1), vessel_center_position(2), vessel_center_position(3));
        normalized_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_min_edge_int(tmp_th_idx,tmp_r_idx)/abs_center_int(tmp_th_idx,tmp_r_idx);
    end
    fprintf('Computing vessels of radius %f. Elapsed time is %f seconds\n', vessel_radius, toc(tmp_tic));
end
%% Save result
edge_int_str = struct;
edge_int_str.psf_FWHM_um = psf_FWHM_um;
edge_int_str.vsl_r_list = vsl_r_list;
edge_int_str.vsl_theta_list = vsl_theta_list;
edge_int_str.psf_sigma_um = psf_sigma;
edge_int_str.simulation_grid_size = dx;
edge_int_str.abs_min_edge_int = abs_min_edge_int;
edge_int_str.abs_center_int = abs_center_int;
edge_int_str.normalized_min_edge_int = normalized_min_edge_int;
edge_int_str.normalized_radial_int_dist = n_radial_int_dist;
edge_int_str.normalized_vertical_radial_int_dist= n_vertical_radial_int_dist;

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