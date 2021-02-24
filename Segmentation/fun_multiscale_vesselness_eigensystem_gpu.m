function [vesselness, eig1, eig2, eig3, A11, A12, A13, ...
    A22, A23, A33] = fun_multiscale_vesselness_eigensystem_gpu(inputArray, vessel_parameters)
% fun_multiscale_vesselness_eigensystem_gpu computes the vesselness of each voxel on
% multiple scale. 
% Input: 
%     inputArray: 3D array
%     vessel_parameters: struct, with field: sigma_list, alpha, beta,
%     gamma, c
% 
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
image_size_padded = size(inputArray) + 2;
vesselness = zeros(image_size_padded);
A11 = zeros(image_size_padded);
A12 = zeros(image_size_padded);
A13 = zeros(image_size_padded);
A22 = zeros(image_size_padded);
A23 = zeros(image_size_padded);
A33 = zeros(image_size_padded);
eig1 = zeros(image_size_padded);
eig2 = zeros(image_size_padded);
eig3 = zeros(image_size_padded);
% Compute the gaussian smoothed image first. 
I_smoothed = cell(numel(vessel_parameters.sigma_list,1));
tmpI = single(gpuArray(inputArray));
for sigma_idx = 1 : numel(vessel_parameters.sigma_list)
        I_smoothed{sigma_idx} = padarray(gather(imgaussfilt3(tmpI, vessel_parameters.sigma_list(sigma_idx), 'FilterDomain', 'spatial')), [1,1,1], 'both', 'replicate');
end
clear tmpI
% Send data to GPU and keep the highest response values. 
for starting_sec_idx = 1 : num_processing_sections : image_size_padded(3)
    tmp_processing_sections = starting_sec_idx : min(image_size_padded(3), starting_sec_idx + num_processing_sections - 1);
    tmp_num_section = length(tmp_processing_sections);
    gpuA11 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuA12 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuA13 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuA22 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuA23 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuA33 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuEig1 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuEig2 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuEig3 = zeros([image_size_padded(1:2), tmp_num_section], 'gpuArray');
    gpuVesselness = zeros([image_size_padded(1:2), tmp_num_section]);
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
        tmp_larger = tmp_vesselness > gpuVesselness;
        gpuEig1 = gpuEig1 .* ~tmp_larger + tmpEig1 .* tmp_larger;
        gpuEig2 = gpuEig2 .* ~tmp_larger + tmpEig1 .* tmp_larger;
        gpuEig3 = gpuEig3 .* ~tmp_larger + tmpEig1 .* tmp_larger;
        gpuA11 = gpuA11 .* ~tmp_larger + tmpA11 .* tmp_larger;
        gpuA12 = gpuA12 .* ~tmp_larger + tmpA12 .* tmp_larger;
        gpuA13 = gpuA13 .* ~tmp_larger + tmpA13 .* tmp_larger;
        gpuA22 = gpuA22 .* ~tmp_larger + tmpA22 .* tmp_larger;
        gpuA23 = gpuA23 .* ~tmp_larger + tmpA23 .* tmp_larger;
        gpuA33 = gpuA33 .* ~tmp_larger + tmpA33 .* tmp_larger;
        gpuVesselness = gpuVesselness .* ~tmp_larger + tmp_vesselness .* tmp_larger;
    end
    A11(:,:,tmp_processing_sections) = gather(gpuA11);
    A12(:,:,tmp_processing_sections) = gather(gpuA12);
    A13(:,:,tmp_processing_sections) = gather(gpuA13);
    A22(:,:,tmp_processing_sections) = gather(gpuA22);
    A23(:,:,tmp_processing_sections) = gather(gpuA23);
    A33(:,:,tmp_processing_sections) = gather(gpuA33);
    eig1(:,:,tmp_processing_sections) = gather(gpuEig1);
    eig2(:,:,tmp_processing_sections) = gather(gpuEig2);
    eig3(:,:,tmp_processing_sections) = gather(gpuEig3);
    vesselness(:,:,tmp_processing_sections) = gather(gpuVesselness);
end
end
