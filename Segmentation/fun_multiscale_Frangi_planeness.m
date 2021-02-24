function planeness = fun_multiscale_Frangi_planeness(inputArray, vessel_parameters)
% fun_multiscale_vesselness_gpu computes the vesselness of each voxel on
% multiple scale. 
% Input: 
%     inputArray: 3D array
%     vessel_parameters: struct, with field: DoG_scale_list, alpha, beta,
%     gamma, c
% 
% Set default values & Convert the input array to single if not float. 
if isa(inputArray, 'gpuArray')
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
    vessel_parameters.DoG_scale_list = [1, 2, 4, 8, 16];
end

if ~isfield(vessel_parameters, 'alpha')
    vessel_parameters.alpha = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'beta')
    vessel_parameters.beta = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'gamma')
    input_max = max(inputArray(:));
    input_min = min(inputArray(:));
    est_bg_eig_sqr_sum = 0.05 * (input_max - input_min);
    vessel_parameters.gamma = 2 * (est_bg_eig_sqr_sum) ^2;
end



image_size = size(inputArray);
planeness = zeros(image_size, 'like', inputArray);
for sigma_idx = 1 : numel(vessel_parameters.DoG_scale_list)
    if vessel_parameters.DoG_scale_list(sigma_idx) > 0
        I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale_list(sigma_idx), 'FilterDomain', 'spatial');
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
    tmp_planeness = exp(- ((tmpEig2./tmpEig3).^2 ) ./ vessel_parameters.alpha);
    tmp_planeness = tmp_blob_measure .* tmp_planeness .* tmp_signal_measure;
    tmp_planeness(isnan(tmp_planeness)) = 0;
    % Remove dark plane
    tmp_planeness = tmp_planeness - tmp_planeness .* ( tmpEig3  >= 0 );
    planeness = max(planeness, tmp_planeness);
end
end


