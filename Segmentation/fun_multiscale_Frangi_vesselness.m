function [tmp_vesselness]= fun_multiscale_Frangi_vesselness(inputArray, vessel_parameters)
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
    vessel_parameters.DoG_scale_list = [1, 2, 4, 8, 16];
end
if ~isfield(vessel_parameters, 'alpha')
    vessel_parameters.alpha = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'beta')
    vessel_parameters.beta = 2 * (0.5) ^2;
end
if ~isfield(vessel_parameters, 'input_normalized')
    vessel_parameters.input_normalized = false;
end
% vessel_parameters.DoG_scale_list = sort(vessel_parameters.DoG_scale_list, 'ascend');
% if nargout > 1
%     record_edge_Q = true;
% else
%     record_edge_Q = false;
% end

image_size = size(inputArray);
vesselness = zeros(image_size, 'like', inputArray);
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

%% Run on CPU - Assume the noise level to be 
%     [tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(double(I_smoothed), vessel_parameters.DoG_scale_list(sigma_idx),1);
%     [tmpEig1, tmpEig2, tmpEig3] = fun_vectorized_3x3_real_sym_mat_eigenvalues(tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);
% %     [tmpEig1, tmpEig2, tmpEig3] = arrayfun(@fun_gpu_3x3_real_sym_mat_eigenvalues,tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33, 0);
%     % deviation from background: 1 means signal
%     tmp_signal_measure = tmpEig1.^2 + tmpEig2.^2 + tmpEig3.^2;
%     % Since most of the image are background, use the average value of the
%     % sum of the squared eigenvalues to estiamte the background
%     tmp_signal_measure = (1 - exp(- (tmp_signal_measure)./mean(tmp_signal_measure(:))));
%     fun_gpu_eigenvalues_to_Frangi
%     tmp_vesselness =  (1 - exp(- ((tmpEig2./tmpEig3).^2 ) ./ vessel_parameters.alpha)) ... % deviation from plane: 1 means plane
%         .* exp(-( tmpEig1.^2 ./abs(tmpEig2 .* tmpEig3) ) ./ vessel_parameters.beta) ... % deviation from blob: 1 means eig1 << eig2 * eig3 -> plane or rod
%         .* tmp_signal_measure;
%     tmp_vesselness( isnan(tmp_vesselness) ) = 0;
%     tmp_vesselness = tmp_vesselness - tmp_vesselness .* (tmpEig2>=0 | tmpEig3 >=0);
%% Completely run on GPU
    if ~vessel_parameters.input_normalized
        I_smoothed = rescale(I_smoothed);
    end
    [tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(I_smoothed, vessel_parameters.DoG_scale_list(sigma_idx),1);
%     if record_edge_Q
%         varargout{1} = gather(sqrt(tmpA11 .^ 2 + tmpA22 .^ 2 + tmpA33 .^ 2));
%         record_edge_Q = false;
%     end
    tmp_vesselness = arrayfun(@fun_gpu_Frangi_vesselness, tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);
%%
    vesselness = max(vesselness, tmp_vesselness);
end
end


