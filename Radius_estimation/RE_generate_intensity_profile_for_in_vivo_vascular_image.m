clear;clc
set_env;
%% Set simulation parameters
DataManager = FileManager;
gpuDevice(2);
vessel_radius_list = [1: 0.1 : 3, 3.2 : 0.2 : 5, 5.5 : 0.5 : 10, 11 : 1 : 15];
z_component_list = 0 : 0.05 : 1;
theta_list = asin(z_component_list);

psf_sigma_xy_um = 0.3750;
psf_sigma_z_um = 2.40;

psf_FWHM_xy_um = psf_sigma_xy_um * 2 * sqrt(2 * log(2));
psf_FWHM_z_um = psf_sigma_z_um * 2 * sqrt(2 * log(2));

bc_layer_thickness_list = [0 : 0.1 : 2, 2.2 : 0.2 : 3];
rbc_int_n_list = [0 : 0.1 : 1];

num_bc_t = numel(bc_layer_thickness_list);
num_rbc_n_int = numel(rbc_int_n_list);

dx_um = 0.1;
%%
invivo_vsl_im_str = struct;
invivo_vsl_im_str.rbc_int_n_list = rbc_int_n_list;
invivo_vsl_im_str.vsl_bc_lyr_thkn_um_list = bc_layer_thickness_list;
invivo_vsl_im_str.vsl_rds_list_um = vessel_radius_list;
invivo_vsl_im_str.vsl_ori_z_comp = z_component_list;
invivo_vsl_im_str.vsl_ele_agl_rad = theta_list;
invivo_vsl_im_str.dx_um = dx_um;
invivo_vsl_im_str.edge_int_str = cell(num_rbc_n_int, num_bc_t);

invivo_vsl_im_str.psf_sigma_xy_um = psf_sigma_xy_um;
invivo_vsl_im_str.psf_sigma_z_um = psf_sigma_z_um;
invivo_vsl_im_str.psf_sigma_um = [psf_sigma_xy_um, psf_sigma_xy_um, psf_sigma_z_um];
invivo_vsl_im_str.psf_FWHM = [psf_FWHM_xy_um, psf_FWHM_xy_um, psf_FWHM_z_um];
invivo_vsl_im_str.psf_FWHM_xy_um = psf_FWHM_xy_um;
invivo_vsl_im_str.psf_FWHM_z_um = psf_FWHM_z_um;
%% Computation
for iter_bc_t = 1 : num_bc_t
    tmp_tic = tic;
    tmp_bc_layer_thickness = bc_layer_thickness_list(iter_bc_t);
    %% Calculate parameters
    tmp_fp = fullfile(DataManager.Scratch_Folder_Path, 'tmp', 'in_vivo_image_simulation', sprintf('psf_data_%d_%d_%d.mat', 1, 1, iter_bc_t));
    if isfile(tmp_fp)
        edge_intensity = load(tmp_fp);
    else
        edge_intensity = fun_radius_estimation_compute_in_vivo_vessel_edge_int_given_psf(invivo_vsl_im_str.psf_sigma_um, ...
            vessel_radius_list, theta_list, tmp_bc_layer_thickness, rbc_int_n_list, invivo_vsl_im_str.dx_um, true);
    end
    % Split the simulation result
    tmp_edge_int_str_cell = cell(num_rbc_n_int, 1);
    for iter_int = 1 : num_rbc_n_int
       tmp_str = edge_intensity;
       tmp_str.vsl_rbc_int_n = tmp_str.vsl_rbc_int_n(iter_int);
       tmp_str.abs_min_edge_int = squeeze(tmp_str.abs_min_edge_int(iter_int, :, :));
       tmp_str.abs_center_int = squeeze(tmp_str.abs_center_int(iter_int, :, :));
       tmp_str.normalized_min_edge_int = squeeze(tmp_str.normalized_min_edge_int(iter_int, :, :));
       tmp_str.normalized_radial_int_dist = squeeze(tmp_str.normalized_radial_int_dist(iter_int, :, :));
       tmp_str.normalized_vertical_radial_int_dist = squeeze(tmp_str.normalized_vertical_radial_int_dist(iter_int, :, :));
       tmp_str.max_normalized_min_edge_int = squeeze(tmp_str.max_normalized_min_edge_int(iter_int, :, :));
       tmp_str.max_normalized_radial_int_dist = squeeze(tmp_str.max_normalized_radial_int_dist(iter_int, :, :));
       tmp_str.max_normalized_vertical_radial_int_dist = squeeze(tmp_str.max_normalized_vertical_radial_int_dist(iter_int, :, :));        
       
       tmp_int = struct;
       tmp_z_comp = sin(theta_list);
       [tmp_int.z_comp, tmp_int.radius] = ndgrid(tmp_z_comp, vessel_radius_list);
       tmp_int.variable_name = {'z_component', 'radius'};
       tmp_int.Interpolation_method = 'linear';
       tmp_int.Extrapolation_method = 'nearest';
       tmp_int.n_min_edge_int = griddedInterpolant(tmp_int.z_comp, tmp_int.radius,...
           tmp_str.normalized_min_edge_int, tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
       tmp_str.n_min_edge_int_interpolation = tmp_int;
       
       tmp_int.n_min_edge_int = griddedInterpolant(tmp_int.z_comp, tmp_int.radius, tmp_str.max_normalized_min_edge_int, ...
           tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
       tmp_str.max_n_min_edge_int_interpolation = tmp_int;
       
       tmp_edge_int_str_cell{iter_int} = tmp_str;
    end
    
    invivo_vsl_im_str.edge_int_str(:, iter_bc_t) = tmp_edge_int_str_cell;
    tmp_folder = fileparts(tmp_fp);
    if ~isfolder(tmp_folder)
        mkdir(tmp_folder);
    end
    save(tmp_fp, '-struct', 'edge_intensity');
    fprintf('Finish the simulation for boundary layer thickness %f um. Elapsed time is %f seconds.\n', ...
        tmp_bc_layer_thickness, toc(tmp_tic));
end
fprintf('Finish in vivo vessel image simulation.\n');
psf_file_fp = fullfile(DataManager.fp_metadata_file('Vessel_radius_calibration', 'In_vivo_image_simulation', 'In_vivo_image_simulation_v2'));
invivo_vsl_im_str.filepath = psf_file_fp;
DataManager.write_data(invivo_vsl_im_str.filepath, invivo_vsl_im_str);

