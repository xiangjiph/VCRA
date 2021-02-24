%% Check 1D case
system_length = 102;
diff_1D_array = false(system_length, 1);
diff_1D_array(1) = true;
diff_1D_array(end) = true;

oxygen_metabolism_rate = 90e-6; %mol/(L sec)
cap_pO2 = 40;
tissue_pO2_0 = 20;
downsample_ratio = 1;

eta = (oxygen_metabolism_rate) * (1/(diffCoeffO2 * 1e12)) * (1/alphaO2); %mmHg/um^2
eta_ds = eta * (downsample_ratio ^2);

target_array_size = round(size(diff_1D_array) ./ downsample_ratio);
vis_recon_rz = imresize(uint8(diff_1D_array), target_array_size, 'Method', 'nearest') > 0;

vis_recon_rz_dt = bwdist(vis_recon_rz) .* downsample_ratio;
vessel_ind = find(vis_recon_rz);
simu_array_size = size(vis_recon_rz);
% Initialization 
simu_array = ones(simu_array_size, 'double') * tissue_pO2_0;
% simu_array = oxygen_metabolism_rate / ( 2 * alphaO2 * diffCoeffO2 * 1e12 ) .* ...
%     ((vis_recon_rz_dt - system_length/2) .^ 2 - system_length^2/4) + cap_pO2 + 1;

simu_array(vessel_ind) = cap_pO2;
simu_array = gpuArray(simu_array);
tic
num_max_iter = 5000;
target_accuracy = 1e-8;
cum_err = zeros(num_max_iter, 1, 'gpuArray');
max_err = 1;
iter = 0;
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
while iter < num_max_iter && max_err > target_accuracy
    iter = iter + 1;
    % Update non-boundary voxels
    tmp_old_value = simu_array(2 : (end - 1));
    tmp_new_vlaue = - eta_ds/2 + ...
       (simu_array(1 : (end - 2)) + simu_array(3 : (end))) / 2;
    simu_array(2 : (end - 1)) = tmp_new_vlaue;
    max_err = max(abs(tmp_new_vlaue - tmp_old_value), [], 'all');
    cum_err(iter) = max_err;
    fprintf('Finish computing iteration %d. Maximum update is %f.\n', ...
        iter, max_err);
    %% Online visualization
    tmp_x_data = 1 : (system_length - 2);
    theo_y = oxygen_metabolism_rate / ( 2 * alphaO2 * diffCoeffO2 * 1e12 ) .* ...
        ((tmp_x_data - (system_length - 1) /2) .^ 2 - (system_length - 1)^2/4) + cap_pO2;
    plt_krogh = plot(ax_hdl, tmp_x_data, theo_y, 'LineWidth', 2, 'LineStyle', ':');
    hold(ax_hdl, 'on');
    tmp_y_data = simu_array(2 : (end-1));
    plt_y_med_hdl = plot(ax_hdl, tmp_x_data, tmp_y_data, 'LineWidth', 2);    
    hold(ax_hdl, 'off');
    pause(0.01);
end
toc
simu_array = gather(simu_array);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
tmp_y_data = simu_array(2 : (end-1));
tmp_x_data = 1 : (system_length - 2);
plt_y_med_hdl = plot(ax_hdl, tmp_x_data, tmp_y_data, 'LineWidth', 2);

theo_y = oxygen_metabolism_rate / ( 2 * alphaO2 * diffCoeffO2 * 1e12 ) .* ...
    ((tmp_x_data - (system_length - 1) /2) .^ 2 - (system_length - 1)^2/4) + cap_pO2;
hold(ax_hdl, 'on');
plt_krogh = plot(ax_hdl, tmp_x_data, theo_y, 'LineWidth', 2, 'LineStyle', ':');
%% 3D sudo-1D case
system_length = 62;
vis_recon = false(system_length, system_length, system_length);
vis_recon(1, :, :) = true;
vis_recon(system_length, :, :) = true;
%% Solve Poisson equation using central finite difference
oxygen_metabolism_rate = 90e-6; %mol/(L sec)
cap_pO2 = 40;
tissue_pO2_0 = 0;
downsample_ratio = 2;
eta_ds = eta * (downsample_ratio ^2);

target_array_size = round(size(vis_recon) ./ downsample_ratio);
vis_recon_rz = imresize3(uint8(vis_recon), target_array_size) > 0;
vis_recon_rz_dt = bwdist(vis_recon_rz) .* downsample_ratio;
vessel_ind = find(vis_recon_rz);
simu_array_size = size(vis_recon_rz);
% Initialization 
simu_array = ones(simu_array_size, 'double') * tissue_pO2_0;

% simu_array = oxygen_metabolism_rate / ( 2 * alphaO2 * diffCoeffO2 * 1e12 ) .* ...
%     ((vis_recon_rz_dt - (system_length - 1) /2) .^ 2 - (system_length - 1)^2/4) + cap_pO2 + 1;

simu_array(vessel_ind) = cap_pO2;
simu_array = gpuArray(simu_array);
tic
num_max_iter = 10000;
target_accuracy = 1e-16;
cum_err = zeros(num_max_iter, 1, 'gpuArray');
max_err = 1;
iter = 0;
while iter < num_max_iter && max_err > target_accuracy
    iter = iter + 1;
    tmp_array_0 = simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1));
    tmp_array_new = - eta_ds/6 + ...
       (simu_array(1 : (end - 2), 2 : (end - 1), 2 : (end - 1)) + ...
        simu_array(3 : (end    ), 2 : (end - 1), 2 : (end - 1)) + ...
        simu_array(2 : (end - 1), 1 : (end - 2), 2 : (end - 1)) + ...
        simu_array(2 : (end - 1), 3 : (end    ), 2 : (end - 1)) + ...
        simu_array(2 : (end - 1), 2 : (end - 1), 1 : (end - 2)) + ...
        simu_array(2 : (end - 1), 2 : (end - 1), 3 : (end    ))) / 6;
    % Update non-boundary voxels
    simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1)) = tmp_array_new;
    % Update boundary voxels
    simu_array(1, :, :) = simu_array(2, :, :);
    simu_array(end, :, :) = simu_array(end - 1, :, :);
    simu_array(:, 1, :) = simu_array(:, 2, :);
    simu_array(:, end, :) = simu_array(:, end - 1, :);
    simu_array(:, :, 1) = simu_array(:, :, 2);
    simu_array(:, :, end) = simu_array(:, :, end - 1);
    % Fix the oxygen partial pressure on the boundary
    simu_array(vessel_ind) = cap_pO2;
    % Compute residual error    
    max_err = max(abs(tmp_array_0 - tmp_array_new), [], 'all');
    cum_err(iter) = max_err;
    fprintf('Finish computing iteration %d. Maximum update is %f.\n', ...
        iter, max_err);
end
toc
simu_array = gather(simu_array);
%% Visualize plane 
vis_sec = 30;
fig_hdl = figure;
ax_hdl_1 = subplot(1,3,1);
imagesc(ax_hdl_1, simu_array(:, :, vis_sec));
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_2 = subplot(1,3,2);
imagesc(ax_hdl_2, -vis_recon_rz_dt(:, :, vis_sec));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_3 = subplot(1,3,3);
imagesc(ax_hdl_3, simu_array(:, :, vis_sec) - simu_array(:, :, 1));
ax_hdl_3.DataAspectRatio = [1,1,1];
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
system_length = 62;
% system_length = system_length / downsample_ratio;
sample_system_length = system_length / downsample_ratio;
sample_ind = sub2ind(target_array_size, (1:sample_system_length)',...
    round(sample_system_length/2 * ones(sample_system_length, 1)), ...
    round(sample_system_length/2 * ones(sample_system_length, 1)));
tmp_x_data = vis_recon_rz_dt(sample_ind);
tmp_y_data = simu_array(sample_ind);

plt_y_med_hdl = plot(ax_hdl, tmp_x_data, tmp_y_data, 'LineWidth', 2);

theo_y = oxygen_metabolism_rate / ( 2 * alphaO2 * diffCoeffO2 * 1e12 ) .* ...
    ((tmp_x_data - (system_length - 2) /2) .^ 2 - (system_length - 2)^2/4) + cap_pO2;
hold(ax_hdl, 'on');
plt_krogh = plot(ax_hdl, tmp_x_data, theo_y, 'LineWidth', 2, 'LineStyle', ':');