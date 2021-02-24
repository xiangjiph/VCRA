function rod_response = fun_max_rod_response(image_array, rod_filter_fun_handle, rod_length, rod_radius, ang_spacing,useGPUQ)
% fun_max_rod_response computes the maximum rod response of each voxel in
% the 3D image array using FFT-based convlution. 
%
%
%
if nargin < 6
    warning('useGPUQ = false by default.')
    useGPUQ = false;
end

theta_range = [0, pi];
phi_range = [0, pi];
block_size = size(image_array);
%% Generate rod filters library
filter_size = rod_length;
num_theta = round(diff(theta_range)/ang_spacing);
num_phi = round(diff(phi_range)/ang_spacing);
% theta = pi/2 is special. 
theta_half_pi = false;
num_angle = 0;
rod_filter_orientation = zeros(2, num_theta * num_phi);
for theta_idx = 1 : num_theta
    theta = theta_range(1) + ang_spacing * (theta_idx - 1);
    for phi_idx = 1 : num_phi
        phi = phi_range(1) + ang_spacing * (phi_idx - 1);
        % Remove duplicate Alternative way is to rotate the starting
        % angle 
        if theta ~= pi/2
            rod_filter_orientation(:,num_angle + 1) = [theta, phi];
            num_angle = num_angle + 1;
        elseif ~theta_half_pi
            rod_filter_orientation(:,num_angle + 1) = [theta, phi];
            num_angle = num_angle + 1;
            theta_half_pi = true;
        end
    end
end

rod_filter_orientation = rod_filter_orientation(:, 1:num_angle);
rod_filter_lib = zeros(filter_size, filter_size, filter_size, num_angle);
for angle_idx = 1 : num_angle
    rod_filter_lib(:,:,:,angle_idx) = rod_filter_fun_handle(rod_length, rod_radius, rod_filter_orientation(1, angle_idx), rod_filter_orientation(2, angle_idx));
end
%% Compute maximum rod response of each voxel
% What we actually want to compute is the correlation of the local pattern
% with the rod-shaped filter of different orientation. The implementation
% of correlation in FFT has a subtile indexing problem and the current
% implementation flip the idx in all dimension. For convienence, here we
% just compute the convlution. 
num_rod_filters = length(rod_filter_orientation);
if useGPUQ
    image_array = gpuArray(image_array);
    rod_filter_lib = gpuArray(rod_filter_lib);
    rod_response = zeros(block_size, 'single','gpuArray');
end
disp('Max rod response');
convFFT_best_size = convFFT_optimized_size(block_size + rod_length - 1);
for rodIdx = 1 : num_rod_filters
    rod_response = max(rod_response, convFFT(image_array, rod_filter_lib(:,:,:,rodIdx), 'same', convFFT_best_size));
end

%% Boundary response correction 
% num_element_coverred = convFFT_v2(num_element_coverred, ones([filter_size,filter_size,filter_size], 'single'), 'same', convFFT_best_size);
% rod_response = rod_response./(num_element_coverred./(filter_size^3)); 


end
