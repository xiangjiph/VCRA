function vesselness= fun_Frangi_vesselness(inputArray, vessel_parameters)
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
if ~isfield(vessel_parameters, 'DoG_scale')
    vessel_parameters.DoG_scale = 1;
end
if ~isfield(vessel_parameters, 'alpha')
    vessel_parameters.alpha = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'beta')
    vessel_parameters.beta = 2 * (0.5) ^2;
end

if vessel_parameters.DoG_scale_list > 0
    if input_gpuArrayQ 
        I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale, 'FilterDomain', 'spatial');
    else
        I_smoothed = imgaussfilt3(inputArray, vessel_parameters.DoG_scale);
    end   
else
    I_smoothed = inputArray;
end
[tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(double(I_smoothed), vessel_parameters.DoG_scale,1);
[tmpEig1, tmpEig2, tmpEig3] =fun_vectorized_3x3_real_sym_mat_eigenvalues(tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);

% deviation from background: 1 means signal
tmp_signal_measure = tmpEig1.^2 + tmpEig2.^2 + tmpEig3.^2;
% Since most of the image are background, use the average value of the
% sum of the squared eigenvalues to estiamte the background
tmp_signal_measure = (1 - exp(- (tmp_signal_measure)./mean(tmp_signal_measure(:))));

vesselness =  (1 - exp(- ((tmpEig2./tmpEig3).^2 ) ./ vessel_parameters.alpha)) ... % deviation from plane: 1 means plane
    .* exp(-( tmpEig1.^2 ./abs(tmpEig2 .* tmpEig3) ) ./ vessel_parameters.beta) ... % deviation from blob: 1 means eig1 << eig2 * eig3 -> plane or rod
    .* tmp_signal_measure;
vesselness = vesselness - vesselness .* (tmpEig2>=0 | tmpEig3 >=0);
vesselness( isnan(vesselness) ) = 0;
end


