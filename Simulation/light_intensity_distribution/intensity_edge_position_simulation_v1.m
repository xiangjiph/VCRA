% %% One dimension:
% % One dimension is a boring case. Analytical calculation shows that the
% % edge is always 1/2 of the center intensity.
% dx = 0.01; % In micron
% x_range = [0, 20];
% num_x = round(diff(x_range)/dx) + 1;
% vessel_radius = 3; % In micron
% vessel_center_position = 10;
% vessel_edge = [vessel_center_position - vessel_radius, vessel_center_position + vessel_radius];
% l_x = x_range(1) : dx : x_range(2);
% vessel_object = zeros(1, num_x);
% vessel_mask = (l_x >= vessel_edge(1)) & ...
%     (l_x <= vessel_edge(2));
% vessel_edge_idx = ceil(vessel_edge/dx) + 1;
% vessel_object(vessel_mask) = 1;
% psf_sigma = 0.25;
% psf_sigma_in_voxel = psf_sigma/dx;
% psf = fun_gaussian_kernel(psf_sigma_in_voxel);
% vessel_image = conv(vessel_object, psf, 'same');
% vessel_image = vessel_image./max(vessel_image(:));
% vessel_edge_int_normalized = vessel_image(vessel_edge_idx);
%
% figure;
% plot(l_x, vessel_object, ':')
% hold on
% plot(l_x, vessel_image)
% title(sprintf('Normalized edge intensity = %f', vessel_edge_int_normalized(1)));
%% Two dimension:
% dx = 0.005; % In micron
% x_range = [0, 16];
% 
% num_x = round(diff(x_range)/dx) + 1;
% [l_x, l_y] = meshgrid( x_range(1) : dx : x_range(2), x_range(1) : dx : x_range(2));
% mask_size = size(l_x);
% vessel_radius = 7; % In micron
% vessel_center_position = ones(1,3) .* x_range(2)/2;
% 
% vessel_center_position_idx = find(l_x == vessel_center_position(1) & l_y == vessel_center_position(2));
% [vc_sub1, vc_sub2] = ind2sub(mask_size, vessel_center_position_idx);
% 
% vessel_object = zeros(num_x, num_x);
% vessel_mask = ((l_x - vessel_center_position(1)).^2 + (l_y - vessel_center_position(2)).^2) < vessel_radius^2;
% vessel_object(vessel_mask) = 1;
% 
% theta = 0:0.005:(2*pi);
% edge_pos1 = round( (vessel_radius * cos(theta) + vessel_center_position(1))/dx ) + 1;
% edge_pos2 = round( (vessel_radius * sin(theta) + vessel_center_position(2))/dx ) + 1;
% vessel_edge_idx = sub2ind(mask_size, edge_pos1, edge_pos2);
% 
% psf_sigma = [0.5,0.5];
% psf_sigma_in_voxel = psf_sigma/dx;
% psf = fun_gaussian_kernel(psf_sigma_in_voxel);
% vessel_image = conv2(conv2(vessel_object, fun_gaussian_kernel(psf_sigma(1)/dx), 'same'),...
%     fun_gaussian_kernel(psf_sigma(2)/dx)', 'same');
% 
% vessel_image = vessel_image./vessel_image(vessel_center_position_idx);
% vessel_edge_int_normalized = vessel_image(vessel_edge_idx);
% 
% edge_pixel_radius = (l_x(vessel_edge_idx) - vessel_center_position(1)).^2 + (l_y(vessel_edge_idx) - vessel_center_position(2)).^2;
% fprintf('Average edge radius: %f\nStandard deviation: %f\n', mean(edge_pixel_radius), std(edge_pixel_radius));
% 
% % figure;
% subplot(2,3,1)
% imshow(vessel_object);
% title('Vessel object');
% subplot(2,3,2)
% imshow(rescale(psf));
% title('Point spread function');
% subplot(2,3,3)
% imshow(vessel_image);
% hold on
% scatter(edge_pos1, edge_pos2, '.', 'g')
% title('Vessel image');
% subplot(2,3,4:6)
% plot(theta, vessel_edge_int_normalized)
% grid on
% title('Normalized light intensity at the edge');
% xlabel('\theta/rad');
% ylabel('Normalized Intensity');
% fprintf('Edge intensity: %f +/- %f\n', mean(vessel_edge_int_normalized), std(vessel_edge_int_normalized))
%% Three dimension:
% Questions:
% 1. How sensitive is the normalized edge intensity with respect to the
% change of point spread function size and shape?
% 2.
% Generate tube
clc;clear;
vessel_radius_list = 1:0.1:6;
theta_list = (0:1:90) .* (pi/180);
psf_sigma = [0.5, 0.5, 2];
dx = 0.05;
psf_sigma_in_voxel = psf_sigma ./ dx;
psf = fun_gaussian_kernel(psf_sigma_in_voxel);
psf_size = size(psf);
psf_center = ceil(size(psf)/2);
abs_min_edge_int = zeros(numel(theta_list), numel(vessel_radius_list));
abs_center_int = zeros(numel(theta_list), numel(vessel_radius_list));
normalized_min_edge_int = zeros(numel(theta_list), numel(vessel_radius_list));
estimated_radius = zeros(numel(theta_list), numel(vessel_radius_list));
radial_int_dist = cell(numel(theta_list), numel(vessel_radius_list));
%%
% for tmp_r_idx = 1 : numel(vessel_radius_list)
%     for tmp_th_idx = 1 : numel(theta_list)
tmp_th_idx = 1;
tmp_r_idx = 3;
        vessel_radius = vessel_radius_list(tmp_r_idx);
        vessel_length = 20;
        theta = theta_list(tmp_th_idx);
        phi = 0;
        vessel_radius_in_voxel = round(vessel_radius / dx);
        vessel_length_in_voxel = round(vessel_length / dx);
        % Generate three dimensional vessel model
        sigmaX = 1000000 * vessel_length_in_voxel;
        sigmaY = vessel_radius_in_voxel;
        sigmaZ = vessel_radius_in_voxel;
        half_filter_size = (vessel_length_in_voxel - 1)/2;
        filter_size = vessel_length_in_voxel;
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
        [X,Y,Z] = meshgrid(gpuArray(-half_filter_size:half_filter_size));
        % Compute the argument in the exponent. Here we do not normalize the
        % gaussian. Initialize as the position of the points
        vessel_object = single(cat(2, X(:),Y(:),Z(:)));
        clear X Y Z
        % Rotation
        vessel_object = - sum(vessel_object' .* (A * vessel_object'),1);
        vessel_object = reshape(vessel_object,filter_size, filter_size, filter_size);
        vessel_object = exp(vessel_object);
        % Get the volume specified by the parameter.
        % 0.001 account for the value drop in the x direction.
        vessel_object = gather(vessel_object > exp(-1));
        
        block_size = size(vessel_object);
        vessel_center_position = ceil(block_size/2);
        
        vessel_image = gather(convn(convn(convn(gpuArray(single(vessel_object)), fun_gaussian_kernel(psf_sigma_in_voxel(1)), 'same'),...
            fun_gaussian_kernel(psf_sigma_in_voxel(2))', 'same'), ...
            reshape(fun_gaussian_kernel(psf_sigma_in_voxel(3)),1,1,[]), 'same'));
        % Find the edge position of the boundary and the intensity on the edge
        edge_mask = vessel_object & ~imerode(vessel_object, strel('sphere', 1));
        voxel_edge_idx = find(edge_mask);
        [edge_pos1, edge_pos2, edge_pos3] = ind2sub(block_size, voxel_edge_idx);
        % Voxel subscripts on the xy plane that pass the center:
        voxel_on_center_xy_plane_Q = (edge_pos3 == vessel_center_position(3)) & edge_pos1 < block_size(1) - psf_size(1)&...
            edge_pos1 > psf_size(1) +1 & edge_pos2 < block_size(2) - psf_size(2) &...
            edge_pos2 > psf_size(2) +1;
        center_plane_pos1 = edge_pos1(voxel_on_center_xy_plane_Q);
        center_plane_pos2 = edge_pos2(voxel_on_center_xy_plane_Q);
        voxel_ind_on_center_xy_plane = sub2ind(block_size, center_plane_pos1, ...
            center_plane_pos2, ones(numel(center_plane_pos1),1)*vessel_center_position(3));
        voxel_edge_int_on_center_xy_plane = vessel_image(voxel_ind_on_center_xy_plane);
        radial_intensity_distribution = vessel_image(:, vessel_center_position(2), vessel_center_position(3));
        radial_pos = (-half_filter_size: half_filter_size)' .* dx;
        radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, radial_pos, radial_intensity_distribution);
        plot(radial_pos, radial_intensity_distribution)
        abs_min_edge_int(tmp_th_idx,tmp_r_idx) = min(voxel_edge_int_on_center_xy_plane);
        abs_center_int(tmp_th_idx,tmp_r_idx) = vessel_image(vessel_center_position(1), vessel_center_position(2), vessel_center_position(3));
        
        normalized_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_min_edge_int(tmp_th_idx,tmp_r_idx)/abs_center_int(tmp_th_idx,tmp_r_idx);
        vessel_center_plane_mask_dt = bwdist(vessel_image(:,:,vessel_center_position(3)) < abs_min_edge_int(tmp_th_idx,tmp_r_idx));
        estimated_radius(tmp_th_idx,tmp_r_idx) = vessel_center_plane_mask_dt(vessel_center_position(1), vessel_center_position(2)) * dx;
%     end
% end
%% Save result
edge_intensity = struct;
edge_intensity.vessel_radius_list = vessel_radius_list;
edge_intensity.theta_list = theta_list;
edge_intensity.psf_size = psf_sigma;
edge_intensity.simulation_grid_size = dx;
edge_intensity.abs_min_edge_int = abs_min_edge_int;
edge_intensity.abs_center_int = abs_center_int;
edge_intensity.normalized_min_edge_int = normalized_min_edge_int;
edge_intensity.estimated_radius = estimated_radius;
save('./Simulation/edge_intensity_psf05053.mat', '-struct', 'edge_intensity');
%% Visualization
% fprintf('Radius %f\tEstimated radius %f\n', vessel_radius, estimated_radius);
% fprintf('Absolute minimum edge intensity: %f\nAbsolute center intensity: %f\nNormalized edge minimum intensity: %f\n', abs_min_edge_int(tmp_th_idx,tmp_r_idx), abs_center_int(tmp_th_idx,tmp_r_idx), normalized_min_edge_int(tmp_th_idx,tmp_r_idx));
%     figure;
%     subplot(2,3,1)
%     imshowpair(vessel_image(:,:,vessel_center_position(3)), vessel_object(:,:,vessel_center_position(3)));
%     title('Vessel image overlaid with vessel model');
%     subplot(2,3,2)
%     slice(psf, psf_center(1), psf_center(2), psf_center(3))
%     title('Point spread function');
%     daspect([1,1,1]);
%     subplot(2,3,3)
%     vessel_center_plane_mask = false(block_size(1,2));
%     voxel_ind_2D = sub2ind(block_size(1:2), center_plane_pos1, center_plane_pos2);
%     vessel_center_plane_mask(voxel_ind_2D) = true;
%     imshowpair(vessel_image(:,:,vessel_center_position(3)), vessel_center_plane_mask);
%     title(sprintf('Estiamted radius: %f\n', estimated_radius));
%     subplot(2,3,4:6)
%     plot(voxel_edge_int_on_center_xy_plane)
%     title(sprintf('Intensity along the vessel edge on the central xy plane. Minimum:%f\n', normalized_min_edge_int(tmp_th_idx,tmp_r_idx)));

figure
plot(theta_list, normalized_min_edge_int(:,1:10:end));
xlabel('Tilted angle/rad');
ylabel('Normalized edge intensity');
title('Normalized minimum edge intensity');
grid on
tmp_lgd = legend(cellstr(num2str(vessel_radius_list(1:10:end)','%.0f')), 'Location', 'northwest');
title(tmp_lgd, 'Vessel radius/\mum');

%figure
%plot(vessel_radius_list, min_edge_int');
%tmp_lgd = legend(num2str(theta_list, '%.2d\n'));
%title(tmp_lgd, 'Tilted angle/rad');
