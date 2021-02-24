function [vesselness, planeness,blobness ]= fun_multiscale_geometry_measures(inputArray, vessel_parameters)
% fun_multiscale_vesselness_gpu computes the vesselness of each voxel on
% multiple scale. 
% Input: 
%     inputArray: 3D array
%     vessel_parameters: struct, with field: DoG_scale_list, alpha, beta,
%     gamma, c
% 
% Set default values
input_gpuArrayQ = isa(inputArray, 'gpuArray');
if input_gpuArrayQ
    if ~isaUnderlying(inputArray, 'float')
        inputArray = single(inputArray);
    end
elseif ~isa(inputArray, 'float')
    inputArray = single(inputArray);
end

if nargin < 2
    vessel_parameters = struct;
end
if ~isfield(vessel_parameters, 'DoG_scale_list')
    vessel_parameters.DoG_scale_list = [1, 2, 4];
end

if ~isfield(vessel_parameters, 'alpha')
    vessel_parameters.alpha = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'beta')
    vessel_parameters.beta = 2 * (0.5) ^2;
end


image_size = size(inputArray);
vesselness = zeros(image_size, 'like', inputArray);
blobness = zeros(image_size, 'like', inputArray);
planeness = zeros(image_size, 'like', inputArray);
for sigma_idx = 1 : numel(vessel_parameters.DoG_scale_list)
    if vessel_parameters.DoG_scale_list(sigma_idx) > 0
        if input_gpuArrayQ
            I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale_list(sigma_idx), 'FilterDomain', 'spatial');
        else
            I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale_list(sigma_idx));
        end
    else
        I_smoothed = inputArray;
    end
    [tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(double(I_smoothed), vessel_parameters.DoG_scale_list(sigma_idx),1);
    [tmpEig1, tmpEig2, tmpEig3] =fun_vectorized_3x3_real_sym_mat_eigenvalues(tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);
    % deviation from blob: 1 means eig1 << eig2 * eig3 -> plane or rod
    tmp_blob_measure = exp(-( tmpEig1.^2 ./abs(tmpEig2 .* tmpEig3) ) ./ vessel_parameters.beta);
    % deviation from background: 1 means signal
    tmp_signal_measure = tmpEig1.^2 + tmpEig2.^2 + tmpEig3.^2;
    % Since most of the image are background, use the average value of the
    % sum of the squared eigenvalues to estiamte the background
    tmp_est_bg_level = mean(tmp_signal_measure(:));
    tmp_signal_measure = (1 - exp(- (tmp_signal_measure)./tmp_est_bg_level));
%     tmp_signal_measure = (1 - exp(- (tmp_eig_sqr_sum)./vessel_parameters.gamma));
    % deviation from plane: 1 means plane
    tmp_plane_measure = exp(- ((tmpEig2./tmpEig3).^2 ) ./ vessel_parameters.alpha);
    
    tmp_vesselness =  (1 - tmp_plane_measure) .*  tmp_blob_measure .* tmp_signal_measure;
    tmp_vesselness( isnan(tmp_vesselness) ) = 0;
    tmp_vesselness = tmp_vesselness - tmp_vesselness .* (tmpEig2>=0 | tmpEig3 >=0);
    
    tmp_planeness = tmp_blob_measure .* tmp_plane_measure .* tmp_signal_measure;
    tmp_planeness(isnan(tmp_planeness)) = 0;
    % Remove black plane
    tmp_planeness = tmp_planeness - tmp_planeness .* ( tmpEig3  >= 0 );
    
    tmp_blobness = (1 - tmp_blob_measure) .* (1 - tmp_plane_measure) .* tmp_signal_measure;
    tmp_blobness(isnan(tmp_blobness)) = 0; 
    
    vesselness = max(vesselness, tmp_vesselness);
    planeness = max(planeness, tmp_planeness);
    blobness = max(blobness, tmp_blobness);
end
end


