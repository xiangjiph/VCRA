function [th, outputMask] = fun_select_threshold_by_cc_stat(inputArray, options)
% The commented part of the script is for computing the surface area and
% volume for each individual connected components. Currently they are not
% used for determineing the threshold. 
% The current version takes about 0.4sec on GPU. For getting the surface area
% and volume for each connected components, it takes ~ 4 seconds on GPU. 



if nargin < 2
    options = struct;
end

if ~isfield(options, 'threshold_list')
    th_list = 0.45 : 0.05 : 0.65; 
else
    th_list = options.threshold_list;
end

if ~isfield(options, 'visualization')
    options.visualization = true;
end

block_lenght = numel(inputArray);
num_th = length(th_list);
th_min = min(th_list);
th_max = max(th_list);
th_diff = th_list(2) - th_list(1);
input_is_gpuArray = isa(inputArray, 'gpuArray');

num_cc = zeros(num_th, 1);
cc_sf_area = single(zeros(num_th, 1, 'like', inputArray));
cc_volume = single(zeros(num_th, 1, 'like', inputArray));
% tmp_avg_cc_sfA_vol_ratio = zeros(num_th, 1);
for tmp_idx = 1 : num_th
    array_mask = inputArray > th_list(tmp_idx);
    if input_is_gpuArray
        [~, num_cc(tmp_idx)] = bwlabeln(gather(array_mask), 26);
    else
        [~, num_cc(tmp_idx)] = bwlabeln(array_mask, 26);
    end

    % Pad array to account for the voxel on the boundary
%     if input_is_gpuArray
%         array_mask = padarray(gpuArray(array_mask), [1,1,1], 0, 'both');
%     else
    array_mask = padarray(array_mask, [1,1,1], 0, 'both');
%     end
%     % Count the volume of each connected components
%     tmp_vol_list = histcounts(array_mask(array_mask>0), num_cc(tmp_idx));
%     % Count the surface area of each connected components
%     % Central difference will give 0 surface area for single isolated voxel
    [tmp_dx, tmp_dy, tmp_dz] = fun_gradient3D(array_mask,[1,2,3], 'intermediate'); 
    tmp_surf_list = nnz(tmp_dx) + nnz(tmp_dy) + nnz(tmp_dz);
%     % The following commented scrpt is for computeing the surface area for
%     % each individual connected components. 
%     tmp_dx = abs(tmp_dx);
%     tmp_dy = abs(tmp_dy);
%     tmp_dz = abs(tmp_dz);
%     tmp_surf_list = histcounts(tmp_dx(tmp_dx>0), tmp_num_cc(tmp_idx)) + ...
%         histcounts(tmp_dy(tmp_dy>0), tmp_num_cc(tmp_idx)) + ...
%         histcounts(tmp_dz(tmp_dz>0), tmp_num_cc(tmp_idx));
%     non_tiny_cc_Q = tmp_vol_list >= 10;
%     tmp_avg_cc_sfA_vol_ratio(tmp_idx) = gather(mean(tmp_surf_list(non_tiny_cc_Q).^3 ./ tmp_vol_list(non_tiny_cc_Q).^2));
    cc_volume(tmp_idx) = nnz(array_mask);
%     cc_volume(tmp_idx) = sum(tmp_vol_list);
    cc_sf_area(tmp_idx) = sum(tmp_surf_list);
    clear tmp_dx tmp_dy tmp_dz array_mask tmp_surf_list
end


% Threshold determined from minimum number of connected components chagne
d_num_cc_abs = abs(diff(num_cc));
d_num_cc_abs_r = d_num_cc_abs./num_cc(1:end-1);
d_th = movmean(th_list, 2, 'Endpoints', 'discard');
th_num_cc = median(d_th(d_num_cc_abs_r == min(d_num_cc_abs_r(:))));
% Threshold determined from local maximin average surface area
avg_surf_area = cc_sf_area./num_cc;
th_max_avg_surface = median(th_list(avg_surf_area  == max(avg_surf_area )));
% Threshold determined from local maximin average volume
avg_vol = cc_volume./num_cc;
th_max_avg_vol = median(th_list(avg_vol== max(avg_vol)));
% Threshold determined form local minimum surface area
surf_vol_ratio_dless = cc_sf_area.^3  ./ cc_volume.^2;
th_min_surf_vol_ratio_dless = median(th_list(surf_vol_ratio_dless == min(surf_vol_ratio_dless)));
% Threshold determined from minimum volume change
d_vol = abs(diff(cc_volume));
th_d_vol = median(d_th(d_vol == min(d_vol(:))));




th_candidate = [th_num_cc, th_max_avg_surface, th_max_avg_vol, th_min_surf_vol_ratio_dless, th_d_vol];
th_between_limitQ = th_candidate > (th_min + th_diff/2 + eps )& th_candidate < (th_max- th_diff/2 - eps);
if ~all(th_between_limitQ)
    warning('Optimized threshold(s) equals the maximum/minimum value of the given threshold list.');
    if nnz(th_between_limitQ)
        disp(th_candidate(th_between_limitQ))
        th = median(th_candidate(th_between_limitQ));
    else
        warning('All the thresholds area equal to the threshold limit');
        th = (th_max + th_min)/2;
    end
else
    th = median(th_candidate);
end

if nargout >1 
    outputMask = inputArray > th;
end

if options.visualization
    figure;
    subplot_num = 5;
    subplot(subplot_num ,1,1)
    plot(d_th,d_num_cc_abs_r);
    title(sprintf('Ratio of connected components change slowest at %.3f',th_num_cc));
    subplot(subplot_num ,1,2)
    plot(th_list, cc_sf_area./num_cc);
    hold on
    plot(th_list, cc_volume./num_cc);
    legend('Average Surface area', 'Average Volume')
    title(sprintf('Average cc surface area and volume maximized at %.3f and %.3f', th_max_avg_surface, th_max_avg_vol));
    subplot(subplot_num ,1,3);
    plot(th_list, cc_volume./block_lenght);
    title('Mask volume ratio');
    % subplot(subplot_num ,1,5)
    % plot(th_list, tmp_avg_cc_sfA_vol_ratio);
    % title('Area^3/Volume^2');
    subplot(subplot_num ,1,4)
    plot(d_th, d_vol);
    title(sprintf('CC volume change rate minimized at %.3f', th_d_vol));
    subplot(subplot_num ,1,5)
    plot(th_list, cc_sf_area.^3 ./ cc_volume.^2);
    set(gca, 'YScale', 'log');
    title(sprintf('(Total Area)^3/ (Total Volume)^2 minimized at %.3f', th_min_surf_vol_ratio_dless));
end

end

