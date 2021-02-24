function block_data = fun_downsample_image_cube(block_data, target_size, use_GPU_Q) %#ok<INUSD>

if nargin < 2
    use_GPU_Q = false;
end
% Fill in the region with 0 pixel value by the median of their neighbor
is_zero_Q = (block_data == 0);
if any(is_zero_Q, 'all')
   is_zero_Q_dilate = imdilate(is_zero_Q, strel('sphere', 1));
   fill_value = median(block_data(is_zero_Q_dilate & ~is_zero_Q));
   block_data(is_zero_Q) = fill_value;
end
block_data = imresize3(block_data, target_size);

% if use_GPU_Q
%     block_data = single(gpuArray(medfilt3(block_data)));
%     block_data = gather(imgaussfilt3(block_data, downsample_smooth_gaussian_sigma, 'FilterDomain', 'spatial'));
%     block_data = output_format(imresize3(block_data, 0.5));
% else
%     block_data = single(medfilt3(block_data));
%     block_data = imgaussfilt3(block_data, downsample_smooth_gaussian_sigma);
%     block_data = output_format(imresize3(block_data, 0.5));
% end




end