function rf = rod_filter(rod_length, radius, theta, phi)
% rod_filter creates a set of rod-shaped filters
% Input: 
%     rod_length: scalar float, lenght of the rod
%     radius: scalar float, radius of the rod(cyliner)
%     theta:  scalar float, rotation angle around the Y axis, starting from the X axis(
%     pi/2 - theta in the sperical coordinate)
%     phi: scalar float, the rotation angle around the Z axis, starting from the X axis
%     (same phi angle in spherical coordiante)
% Output: 
%     rf: rod_length-by-rod_length-by-rod_length double array, rod filter
% Modification 10/01/2018
% 1. Correct the factor of 2 error. Previous implementation gives rod of
% radius sqrt(2) * radius.
sigmaX = 1000*rod_length;
sigmaY = radius;
sigmaZ = radius;

% half_filter_size = round(length/2);
% filter_size = 2 * half_filter_size + 1;
% Modification: 
half_filter_size = (rod_length - 1)/2;
filter_size = rod_length;

% Rotation matrix around Y axis
rotateY = [[cos(theta), 0, sin(theta)];[0,1,0];[-sin(theta), 0, cos(theta)]];
% Rotation matrix around Z axis
rotateZ = [[cos(phi), sin(phi),0];[-sin(phi), cos(phi),0];[0,0,1]];
% Rotate around the Z axis first, followed by the rotation around Y axis. 
rotate3D = rotateY * rotateZ ;
% inverse or the covariance matrix of the 3D gaussian distribution
covMat = [[1/(sigmaX^2), 0, 0]; [0, 1/(sigmaY^2), 0]; [0, 0, 1/(sigmaZ^2)]];
% Rotate the coordinates
A = inv(rotate3D) * covMat * rotate3D;

[X,Y,Z] = meshgrid(-half_filter_size:half_filter_size);
R = cat(2, X(:),Y(:),Z(:));
% Compute the argument in the exponent. Here we do not normalize the
% gaussian
arg = - sum(R' .* (A * R'),1);
arg = reshape(arg,filter_size, filter_size, filter_size);
rf = exp(arg);
% To make sure that the volume of the rod of different orientation are comparialbe, 
% we selece the area that lies within the sphere of radius length/2
filter_sphere_mask = (X.^2 + Y.^2 + Z.^2) < (rod_length/2)^2;
% Get the volume specified by the parameter. 
% 0.001 account for the value drop in the x direction. 
filter_value_mask = rf < exp(-1 + 0.001);
rf(filter_value_mask) = 0;
rf(~filter_sphere_mask) = 0;
rf = rf./sum(rf(:));

end