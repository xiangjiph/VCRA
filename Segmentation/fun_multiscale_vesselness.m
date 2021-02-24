function vesselness = fun_multiscale_vesselness(inputArray, vessel_parameters)
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
    vessel_parameters.DoG_scale_list = [1, 2, 4, 8, 16, 32];
end

if ~isfield(vessel_parameters, 'alpha')
    vessel_parameters.alpha = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'beta')
    vessel_parameters.beta = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'gamma')
    vessel_parameters.gamma = 2 * (5) ^2;
end
if ~isfield(vessel_parameters, 'c')
    vessel_parameters.c =  1e-6;
end
image_size = size(inputArray);
vesselness = zeros(image_size);
for sigma_idx = 1 : numel(vessel_parameters.DoG_scale_list)
    if vessel_parameters.DoG_scale_list(sigma_idx) > 0
        I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale_list(sigma_idx), 'FilterDomain', 'spatial');
    else
        I_smoothed = inputArray;
    end        
    [tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(double(I_smoothed), vessel_parameters.DoG_scale_list(sigma_idx),1);
    [tmpEig1, tmpEig2, tmpEig3] =fun_vectorized_3x3_real_sym_mat_eigenvalues(tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);
    tmp_vesselness = single((1 - exp(- ((tmpEig2./tmpEig3).^2) ./ vessel_parameters.alpha)) ...
        .* exp(-( tmpEig1.^2 ./abs(tmpEig2 .* 2) ) ./ vessel_parameters.beta) ...
        .* (1 - exp(- (tmpEig1.^2 + tmpEig2.^2 + tmpEig3.^2)./vessel_parameters.gamma)) ...
        .* exp(- vessel_parameters.c ./(abs(tmpEig2).*tmpEig3.^2)));
    tmp_vesselness( tmpEig2>=0 | tmpEig3 >=0 ) = 0;
    tmp_vesselness( isnan(tmp_vesselness) ) = 0;
    vesselness = max(vesselness, tmp_vesselness);
end
end


