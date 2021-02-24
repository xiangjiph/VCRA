function [lib, orientation] = rod_filter_library(rod_filter_handle, rod_length, rod_radius, omega_step, theta_range, phi_range, remove_duplicateQ)
% Input:
%   ROD_LENGTH: length of rod, in the number of pixels. 
%   ROD_RADIUS: radius of rod, in the number of pixels. Defined as the value
%       that the gaussian distribution dropped by 1/e
%   OMEGA_STEP: step in the rotation angle, same for both phi and theta. 
%   REMOVE_DUPLICATEQ: When theta = pi/2, any phi gives identical matrix.
%   if REMOVE_DUPLICATEQ = true, remove the duplicate angle. 
% If theta_range and phi_range is not given, use default value [0, pi]

if nargin < 5
    theta_range = [omega_step/2, pi];
    phi_range = [omega_step/2, pi];
    remove_duplicateQ = true;
elseif nargin < 7
    remove_duplicateQ = true;
end

filter_size = rod_length;
num_theta = round(diff(theta_range)/omega_step);
num_phi = round(diff(phi_range)/omega_step);
% theta = pi/2 is special. 
theta_half_pi = false;
num_angle = 0;
orientation = zeros(2, num_theta * num_phi);
for theta_idx = 1 : num_theta
    theta = theta_range(1) + omega_step * (theta_idx - 1);
    for phi_idx = 1 : num_phi
        phi = phi_range(1) + omega_step * (phi_idx - 1);
        if remove_duplicateQ
            if theta ~= pi/2
                orientation(:,num_angle + 1) = [theta, phi];
                num_angle = num_angle + 1;
            elseif ~theta_half_pi
                orientation(:,num_angle + 1) = [theta, phi];
                num_angle = num_angle + 1;
                theta_half_pi = true;
            end
        else
            orientation(:,num_angle + 1) = [theta, phi];
            num_angle = num_angle + 1;
        end
    end
end

orientation = orientation(:, 1:num_angle);
lib = zeros(filter_size, filter_size, filter_size, num_angle);
for angle_idx = 1 : num_angle
    lib(:,:,:,angle_idx) = rod_filter_handle(rod_length, rod_radius, orientation(1, angle_idx), orientation(2, angle_idx));
end


end

