clear;clc
set_env;
%% Set simulation parameters
DataManager = FileManager;
gpuDevice(1);
vessel_radius_list = [0.5: 0.1 : 2.5, 2.75 : 0.25 : 6];
z_component_list = 0 : 0.05 : 1;
theta_list = asin(z_component_list);

psf_sigma_xy_list = sort([0.2 : 0.05 : 0.6, 0.275 : 0.05 : 0.425], 'ascend');
psf_sigma_z_list = sort([1.65 : 0.25 : 4.5, 3.275 : 0.25 : 4.025 ], 'ascend');

psf_FWHM_xy_list = psf_sigma_xy_list * 2 * sqrt(2 * log(2));  
psf_FWHM_z_list = psf_sigma_z_list * 2 * sqrt(2 * log(2));
num_psf_xy = numel(psf_FWHM_xy_list);
num_psf_z = numel(psf_FWHM_z_list);

dx_min_scale_ratio = 0.1;
%%
psf_scan_str = struct;
psf_scan_str.psf_sigma_xy_list = psf_sigma_xy_list;
psf_scan_str.psf_sigma_z_list = psf_sigma_z_list;
psf_scan_str.psf_sigma = cell(num_psf_xy, num_psf_z);
psf_scan_str.psf_FWHM = cell(num_psf_xy, num_psf_z);
psf_scan_str.psf_FWHM_xy_list = psf_FWHM_xy_list;
psf_scan_str.psf_FWHM_z_list = psf_FWHM_z_list;
psf_scan_str.vessel_radius_list_um = vessel_radius_list;
psf_scan_str.vessel_ori_z_comp = z_component_list;
psf_scan_str.vessel_elevation_agl_rad = theta_list;
psf_scan_str.dx_min_scale_ratio = dx_min_scale_ratio;
psf_scan_str.edge_int_str = cell(num_psf_xy, num_psf_z);
%% Computation
for iter_psf_z = num_psf_z : -1 : 1
    for iter_psf_xy = 1 : 1 : num_psf_xy
        tmp_tic = tic;
        psf_FWHM = [ones(1,2) .* psf_FWHM_xy_list(iter_psf_xy), psf_FWHM_z_list(iter_psf_z)];
        psf_sigma = [ones(1,2) .* psf_sigma_xy_list(iter_psf_xy), psf_sigma_z_list(iter_psf_z)];
        %% Calculate parameters
        tmp_fp = fullfile(DataManager.Scratch_Folder_Path, 'tmp', 'psf_fitting_data', sprintf('psf_data_%d_%d.mat', iter_psf_xy, iter_psf_z));
        if isfile(tmp_fp)
            edge_intensity = load(tmp_fp);
        else
            edge_intensity = fun_radius_estimation_compute_vessel_edge_int_given_psf(psf_sigma, ...
                vessel_radius_list, theta_list, psf_scan_str.dx_min_scale_ratio);
        end        
        psf_scan_str.edge_int_str{iter_psf_xy, iter_psf_z} = edge_intensity;
        psf_scan_str.psf_FWHM{iter_psf_xy, iter_psf_z} = psf_FWHM;
        psf_scan_str.psf_sigma{iter_psf_xy, iter_psf_z} = psf_sigma;
        tmp_folder = fileparts(tmp_fp);
        if ~isfolder(tmp_folder)
            mkdir(tmp_folder);
        end
        save(tmp_fp, '-struct', 'edge_intensity');
        fprintf('Finish the simulation for PSF (%f, %f, %f). Elapsed time is %f seconds.\n', ...
            psf_FWHM, toc(tmp_tic));
    end
end
fprintf('Finish PSF parameter swap simulation.\n');
save(psf_scan_str.filepath, '-struct', 'psf_scan_str');

psf_file_fp = fullfile(DataManager.fp_metadata_file('WholeBrain', 'ML_2018_08_15', 'psf_fitting_data'));
psf_scan_str.filepath = psf_file_fp;
save(psf_scan_str.filepath, '-struct', 'psf_scan_str');
DataManager.save_to_server_from_root(psf_scan_str.filepath);

