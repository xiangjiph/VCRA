function pO2_result = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    target_downsample_rate, inhomogeneous_term, mask_itp_method, use_gpu_Q)
% This function solve the 3D Poisson equation with constant inhomogeneous
% term, homogeneous boundary condition on the input geometriy and adiabatic
% boundary condition on the system edge. 
% mask_itp_method = 'nearest';
if nargin < 4
    mask_itp_method = 'linear';
    use_gpu_Q = true;
elseif nargin < 5
    use_gpu_Q = true;
end
assert(mod(log2(target_downsample_rate), 1) == 0, ...
    'Target downsample rate should be {1, 2, 4, 8, 16, 32}');
% Default setting 
downsample_rate_list = [32, 16, 8, 4, 2, 1];
downsample_rate_list = downsample_rate_list(downsample_rate_list >= target_downsample_rate);
% Construct the ROI mask. Test convergence based on the voxel in the ROI
system_size = size(vessel_recon_mask);
roi_mask = fun_downsample_by_block_operation(vessel_recon_mask, @max, [128 128 128], true);
has_empty_block_Q = ~all(roi_mask, 'all');
roi_mask = imresize3(uint8(roi_mask), system_size, 'Method', mask_itp_method) >= 0.5;

if ~has_empty_block_Q
    downsample_rate_list = downsample_rate_list(downsample_rate_list <= 8);
end
num_scale = numel(downsample_rate_list);
%% Multigrid method
mg_solver_tic = tic;
pO2_result = [];
for iter_scale = 1 : num_scale
   tmp_ds_rate = downsample_rate_list(iter_scale);
   if tmp_ds_rate <= 4
       tmp_mask_itp_method = mask_itp_method;
   else
       % When downsampling the image by a facter greater than 4, only
       % nearest neighbor interpolation gives reasonable result, as the
       % typical radius of capillary is only 2 um. 
       tmp_mask_itp_method = 'nearest';
   end
   
   tmp_target_size = round(system_size ./ tmp_ds_rate);
   % After testing several methods, it seems that uint8 + linear at 2 um
   % resolution give the closest DT statistics to the 1 um calculation. 
   if tmp_ds_rate ~= 1
       tmp_recon_rz = imresize3(uint8(vessel_recon_mask), tmp_target_size, 'Method', tmp_mask_itp_method) >= 0.5;
   else
       tmp_recon_rz = vessel_recon_mask;
   end
   if all(tmp_recon_rz, 'all')
       continue;
   end
   if has_empty_block_Q
       if tmp_ds_rate ~= 1
           tmp_roi_mask = imresize3(uint8(roi_mask), tmp_target_size, 'Method', tmp_mask_itp_method) >= 0.5;
       else
           tmp_roi_mask = roi_mask;
       end
   else
       tmp_roi_mask = [];
   end
   tmp_inhomogeneous_term = inhomogeneous_term * tmp_ds_rate ^2;
   if isempty(pO2_result)
       tmp_ini_pO2_n = 0;
   else
       tmp_ini_pO2_n = imresize3(pO2_result.pO2_array, tmp_target_size, 'Method', 'linear');
       tmp_ini_pO2_n(tmp_recon_rz) = 0;
   end 
   pO2_result = fun_simulation_OT_solve_ct_diff_itr_cvg_in_roi(tmp_recon_rz, ...
       tmp_inhomogeneous_term, tmp_ini_pO2_n, tmp_roi_mask, use_gpu_Q);       
end
pO2_result.downsample_rate = target_downsample_rate;
pO2_result.downsample_method = mask_itp_method;
fprintf('Finish solving Poisson equation using multigrid Jacobi method. Elapse time is %f seconds.\n', ...
    toc(mg_solver_tic));
end