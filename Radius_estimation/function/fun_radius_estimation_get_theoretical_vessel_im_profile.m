function result_str = fun_radius_estimation_get_theoretical_vessel_im_profile(vsl_r_pxl, simu_system_half_size, ...
    vessel_evaluation_angle, psf_ker_x, psf_ker_y, psf_ker_z, use_GPU_Q)

if nargin < 7
    use_GPU_Q = false;
end

assert(~isscalar(simu_system_half_size), 'The function has been modified');
%% Parameters
psf_size = [numel(psf_ker_x), numel(psf_ker_y), numel(psf_ker_z)];
% Orientation angle of the vessel 
theta = vessel_evaluation_angle;
phi = 0;
% Generate three dimensional vessel model
sigmaX = 1000000 * vsl_r_pxl;
sigmaY = vsl_r_pxl;
sigmaZ = vsl_r_pxl;

half_x_size = simu_system_half_size(1);
half_y_size = simu_system_half_size(2);
half_z_size = simu_system_half_size(3);

% Rotation matrix around Y axis
rotateY = [[cos(theta), 0, sin(theta)];[0,1,0];[-sin(theta), 0, cos(theta)]];
% Rotation matrix around Z axis
rotateZ = [[cos(phi), sin(phi),0];[-sin(phi), cos(phi),0];[0,0,1]];
% Rotate around the Z axis first, followed by the rotation around Y axis.
rotate3D = rotateY * rotateZ ;
% inverse or the covariance matrix of the 3D gaussian distribution
covMat = [[1/(sigmaX^2), 0, 0]; [0, 1/(sigmaY^2), 0]; [0, 0, 1/(sigmaZ^2)]];
% Rotate the coordinates
A = (rotate3D \ covMat) * rotate3D;
if use_GPU_Q
    [X,Y,Z] = meshgrid(gpuArray(-half_x_size:half_x_size),gpuArray(-half_y_size: half_y_size), gpuArray(-half_z_size:half_z_size));
else
    [X,Y,Z] = meshgrid((-half_x_size:half_x_size), (-half_y_size: half_y_size), (-half_z_size:half_z_size));
end
% Compute the argument in the exponent. Here we do not normalize the
% gaussian. Initialize as the position of the points
vessel_object = cat(1, X(:)',Y(:)',Z(:)');
block_size = size(X);
clear X Y Z
% Rotation
vessel_object = - sum(vessel_object .* (A * vessel_object),1);
vessel_object = reshape(vessel_object,block_size(1), block_size(2), block_size(3));
vessel_object = exp(vessel_object);
% Get the volume specified by the parameter.
% 0.001 account for the value drop in the x direction.
vessel_object = gather(vessel_object > exp(-1));
vessel_image = convn(convn(convn(gpuArray(single(vessel_object)), psf_ker_x, 'same'),...
    psf_ker_y, 'same'), psf_ker_z, 'same');
if isa(vessel_image, 'gpuArray')
    vessel_image = gather(vessel_image);
end

vessel_center_position = ceil(block_size/2);

result_str.vessel_image = vessel_image;
result_str.vessel_mask = vessel_object;
result_str.vessel_center_sub = vessel_center_position;
end