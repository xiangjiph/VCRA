%% About point spread function
% From Phil 2009 Journal of Neuroscience paper: 
% The point spread function at the depth of 500 microns is about 0.6x0.6x4
% micron ( standard deviation?)

% Objective used in MouseLight project
% 40x/1.3 NA oil-immersion objective (Plan-Apochromat 40x/1.3DIC Carl
% Zeiss, Germany)
% Point spread function: % https://elifesciences.org/articles/10566/figures#fig1s1
% FWHM = 0.45 x 0.45 x 1.33 um
% sigma = FWHM/2 sqrt(2ln(2))

% Questions:
% 1. How sensitive is the normalized edge intensity with respect to the
% change of point spread function size and shape? - Compute the edge
% intensity for different shape of point spread function

clear;clc
psf_expand_coeff = 1.6;
psf_FWHM = [0.45 * psf_expand_coeff,0.45 * psf_expand_coeff,1.33 * psf_expand_coeff^2];
vessel_radius_list = [0.5:0.1:0.9, 1:0.25:15];
theta_list = (0:5:90) .* (pi/180);
% psf_sigma = [0.6, 0.6, 4];
psf_sigma = psf_FWHM./(2*sqrt(2*log(2)));
save_filename = sprintf('./Simulation/edge_intensity_%s.mat', num2str(round(psf_sigma*100), '%.3d'));
dx = 0.025;
psf_sigma_in_voxel = round(psf_sigma ./ dx);
psf = fun_gaussian_kernel(psf_sigma_in_voxel);
ker1 = fun_gaussian_kernel(psf_sigma_in_voxel(1));
ker2 = fun_gaussian_kernel(psf_sigma_in_voxel(2))';
ker3 = reshape(fun_gaussian_kernel(psf_sigma_in_voxel(3)),1,1,[]);
psf_size = size(psf);
psf_center = ceil(size(psf)/2);
abs_min_edge_int = zeros(numel(theta_list), numel(vessel_radius_list));
abs_center_int = zeros(numel(theta_list), numel(vessel_radius_list));
normalized_min_edge_int = zeros(numel(theta_list), numel(vessel_radius_list));
estimated_radius = zeros(numel(theta_list), numel(vessel_radius_list));
radial_int_dist = cell(numel(theta_list), numel(vessel_radius_list));
%%
% tmp_th_idx = 1;
% tmp_r_idx = 10;

for tmp_r_idx = numel(vessel_radius_list) : -1 : 1
    tic
    vessel_radius = vessel_radius_list(tmp_r_idx);
    fprintf('Computing vessels of radius %f\n', vessel_radius);
    for tmp_th_idx = 1 : numel(theta_list)
        theta = theta_list(tmp_th_idx);
        phi = 0;
        vessel_radius_in_voxel = round(vessel_radius / dx);
        vessel_length_in_voxel = psf_size(3);
        % Generate three dimensional vessel model
        sigmaX = 1000000 * vessel_length_in_voxel;
        sigmaY = vessel_radius_in_voxel;
        sigmaZ = vessel_radius_in_voxel;
        half_x_size = single(vessel_radius_in_voxel + ceil(psf_size(1)/2)) + 1;
        half_y_size = single(vessel_radius_in_voxel + ceil(psf_size(2)/2)) + 1;
        half_z_size = single(ceil(vessel_length_in_voxel/2)) + 1;

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
        [X,Y,Z] = meshgrid(gpuArray(-half_x_size:half_x_size),gpuArray(-half_y_size: half_y_size), gpuArray(-half_z_size:half_z_size));
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
        vessel_center_position = ceil(block_size/2);
        vessel_image = gather(convn(convn(convn(gpuArray(single(vessel_object)), ker1, 'same'),...
            ker2, 'same'), ker3, 'same'));
        % Find the edge position of the boundary and the intensity on the edge
        edge_mask = vessel_object & ~imerode(vessel_object, strel('sphere', 1));
        voxel_edge_idx = find(edge_mask);
        [edge_pos1, edge_pos2, edge_pos3] = ind2sub(block_size, voxel_edge_idx);
        % Voxel subscripts on the xy plane that pass the center:
        voxel_on_center_xy_plane_Q = (edge_pos3 == vessel_center_position(3)) & edge_pos1 < block_size(1) - psf_size(1)/2&...
            edge_pos1 > psf_size(1)/2 +1 & edge_pos2 < block_size(2) - psf_size(2)/2 &...
            edge_pos2 > psf_size(2)/2 +1;
        center_plane_pos1 = edge_pos1(voxel_on_center_xy_plane_Q);
        center_plane_pos2 = edge_pos2(voxel_on_center_xy_plane_Q);

        voxel_ind_on_center_xy_plane = sub2ind(block_size, center_plane_pos1, ...
            center_plane_pos2, ones(numel(center_plane_pos1),1)*vessel_center_position(3));
        % Edge intensity on the center xy plane
        voxel_edge_int_on_center_xy_plane = vessel_image(voxel_ind_on_center_xy_plane);
        % Intensity distribution along the radius
        radial_intensity_distribution = vessel_image(vessel_center_position(1):end , vessel_center_position(2), vessel_center_position(3));
        abs_lateral_edge_int = radial_intensity_distribution(vessel_radius_in_voxel+1);
        
        radial_pos = (0:(half_x_size))' .* dx;
        radial_int_dist{tmp_th_idx, tmp_r_idx} = cat(2, radial_pos, radial_intensity_distribution);
%         abs_min_edge_int(tmp_th_idx,tmp_r_idx) = min(voxel_edge_int_on_center_xy_plane);
        abs_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_lateral_edge_int;
        abs_center_int(tmp_th_idx,tmp_r_idx) = vessel_image(vessel_center_position(1), vessel_center_position(2), vessel_center_position(3));
        normalized_min_edge_int(tmp_th_idx,tmp_r_idx) = abs_min_edge_int(tmp_th_idx,tmp_r_idx)/abs_center_int(tmp_th_idx,tmp_r_idx);
        % Mask the vessel image on the center xy plane by the minimum
        % edge intensity, compute the distance transform and get the
        % estimated radius
        vessel_center_plane_mask_dt = bwdist(vessel_image(:,:,vessel_center_position(3)) <= abs_min_edge_int(tmp_th_idx,tmp_r_idx));
        estimated_radius(tmp_th_idx,tmp_r_idx) = vessel_center_plane_mask_dt(vessel_center_position(1), vessel_center_position(2)) * dx;
    end
    toc
end

% Save result
edge_intensity = struct;
edge_intensity.vessel_radius_list = vessel_radius_list;
edge_intensity.theta_list = theta_list;
edge_intensity.psf_size = psf_sigma;
edge_intensity.simulation_grid_size = dx;
edge_intensity.abs_min_edge_int = abs_min_edge_int;
edge_intensity.abs_center_int = abs_center_int;
edge_intensity.normalized_min_edge_int = normalized_min_edge_int;
edge_intensity.estimated_radius = estimated_radius;
edge_intensity.radial_int_distribution = radial_int_dist;
save(save_filename, '-struct', 'edge_intensity');
% Visualization
%     figure;
%     subplot(2,3,1)
%     imshowpair(vessel_image(:,:,vessel_center_position(3)), vessel_object(:,:,vessel_center_position(3)));
%     title('Vessel image overlaid with vessel model');
%     subplot(2,3,2)
%     slice(psf, psf_center(1), psf_center(2), psf_center(3))
%     xlabel('x/dx');
%     ylabel('y/dx');
%     zlabel('z/dx');
%     title(sprintf('Point spread function\n[%.1f, %.1f, %.1f]um', psf_sigma(1), psf_sigma(2), psf_sigma(3)));
%     daspect([1,1,1]);
%     subplot(2,3,3)
%     vessel_center_plane_mask = false(block_size(1,2));
%     voxel_ind_2D = sub2ind(block_size(1:2), center_plane_pos1, center_plane_pos2);
%     vessel_center_plane_mask(voxel_ind_2D) = true;
%     imshowpair(vessel_image(:,:,vessel_center_position(3)), vessel_center_plane_mask);
%     title(sprintf('Estiamted radius: %f\n', estimated_radius(tmp_th_idx,tmp_r_idx)));
%     subplot(2,3,4:5)
%     plot(voxel_edge_int_on_center_xy_plane)
%     xlabel('Edge voxel index');
%     ylabel('Normalized edge intensity');
%     title(sprintf('Intensity along the vessel edge on the central xy plane. Minimum:%f\n', abs_lateral_edge_int));
%     subplot(2,3,6)
%     plot(radial_pos, radial_intensity_distribution)
%     grid on
%     xlabel('Radial position/\mum');
%     ylabel('Normalized intensity');
    
%% Plot figure
% figure
% plot(theta_list, normalized_min_edge_int);
% xlabel('Tilted angle/rad');
% ylabel('Normalized edge intensity');
% title('Normalized minimum edge intensity');
% grid on
% tmp_lgd = legend('1.00', '1.25', '1.50', '1.75', '2.00', '3.00', '4.00', '5.00', '6.00');
% title(tmp_lgd, 'Vessel radius/\mum');

figure
plot(vessel_radius_list, normalized_min_edge_int(1:1:end,:)', 'LineWidth', 2);
set(gca, 'LineWidth', 2);
tmp_lgd = legend(cellstr(num2str((0:10:90)', '%.2d\n')));
title(tmp_lgd, 'Tilted angle/degree');
grid on
