function vesselness = fun_multiscale_vesselness_gpu(inputArray, vessel_parameters)
% fun_multiscale_vesselness_gpu computes the vesselness of each voxel on
% multiple scale. 
% Input: 
%     inputArray: 3D array
%     vessel_parameters: struct, with field: sigma_list, alpha, beta,
%     gamma, c
% This function divide the data into small blocks. For array of size 250^3,
% there are no need to do this. 
% Set default values
if nargin < 2
    vessel_parameters = struct;
end
if ~isfield(vessel_parameters, 'sigma_list')
    vessel_parameters.sigma_list = [1, 2, 4, 8, 16, 32];
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
% Number of sections of the image stack sent to GPU for processing in each
% iteration
num_processing_sections = 64;

% Initialization
image_size = size(inputArray);
vesselness = zeros(image_size);
% Compute the gaussian smoothed image first. 
I_smoothed = cell(numel(vessel_parameters.sigma_list,1));
tmpI = single(gpuArray(inputArray));
for sigma_idx = 1 : numel(vessel_parameters.sigma_list)
        I_smoothed{sigma_idx} = gather(imgaussfilt3(tmpI, vessel_parameters.sigma_list(sigma_idx), 'FilterDomain', 'spatial'));
end
clear tmpI
% Send data to GPU and keep the highest response values. 
for starting_sec_idx = 1 : num_processing_sections : image_size(3)
    tmp_processing_sections = starting_sec_idx : min(image_size(3), starting_sec_idx + num_processing_sections - 1);
    tmp_num_section = length(tmp_processing_sections);
    gpuVesselness = zeros([image_size(1:2), tmp_num_section], 'gpuArray');
    for sigma_idx = 1 : numel(vessel_parameters.sigma_list)
        I_sections = gpuArray(I_smoothed{sigma_idx}(:,:, tmp_processing_sections));
        [tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33] = fun_gamma_normalized_hessian_matrix3D(double(I_sections), vessel_parameters.sigma_list(sigma_idx),1);
        [tmpEig1, tmpEig2, tmpEig3] =fun_vectorized_3x3_real_sym_mat_eigenvalues(tmpA11, tmpA12, tmpA13, tmpA22, tmpA23, tmpA33);
        tmp_vesselness = single((1 - exp(- ((tmpEig2./tmpEig3).^2) ./ vessel_parameters.alpha)) ...
            .* exp(-( tmpEig1.^2 ./abs(tmpEig2 .* 2) ) ./ vessel_parameters.beta) ...
            .* (1 - exp(- (tmpEig1.^2 + tmpEig2.^2 + tmpEig3.^2)./vessel_parameters.gamma)) ...
            .* exp(- vessel_parameters.c ./(abs(tmpEig2).*tmpEig3.^2)));
        tmp_vesselness( tmpEig2>=0 | tmpEig3 >=0 ) = 0;
        tmp_vesselness( isnan(tmp_vesselness) ) = 0;
        gpuVesselness = max(gpuVesselness, tmp_vesselness);
    end
    vesselness(:,:,tmp_processing_sections) = gather(gpuVesselness);
end
end


