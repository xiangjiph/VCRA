function [mean_array, std_array] = fun_masked_local_meanNstd(input_array, window_size, mask, out_mask_value,save_to_RAM_Q) 
% fun_masked_local_meanNstd computes the local mean and standard deviation
% for each voxel in the mask. The voxels out of the mask are left out. 
% Input: 
%   input_array: 3 dimensional numerical array. Usually is uint8 or uint16
%   window_size: numerical scalar or 3-by-1 numerical array, specify the
%   size of the window for computing local statistics. 
%   mask: 3 dimensional logical array of the same size as input_array. The
%   voxel out of the mask are set to 0 and not used for computing local
%   statistics. 
%   out_mask_value: {0, 'inf', 'nan'} the value of the out-of-mask voxels
%   save_to_RAM_Q: save the final result to memory ( only functions if
%   input_array is a gpuArray)
% Output:
%   mean_array: 3 dimensional numerical array of the same size as
%   input_array. Mean of the voxel value within both the mask and the
%   window. 
%   std_array: 3 dimensional numerical array of the same size as
%   input_array. Standard deviation of the voxel value within both the mask
%   and the window. *** Normalized by N, instead of (N - 1) ***
if nargin < 4
    save_to_RAM_Q = true;
    out_mask_value = nan;
elseif nargin < 5
    save_to_RAM_Q = true;
end

if isscalar(window_size)
    window_size = ones(3,1) * window_size;
end
% Assume the minimum value of the input_array is 0. Shift the input_array
% by +1 to distinguish between voxel in the mask and out of the mask
% compute_mask = imerode(mask, strel('cube',2));
input_array = single(input_array) + 1;
input_array(~mask) = 0;

ker1 = ones(window_size(1),1, 'like', input_array);
ker2 = ones(1, window_size(2), 'like', input_array);
ker3 = ones(1, 1, window_size(3), 'like', input_array);

% Number of local element
num_elem_array = convn(convn(convn(mask, ker1, 'same'), ker2, 'same'), ker3, 'same');
% Add num_elem_array by 1 to avoid divided-by-zero error
num_elem_array = max(1, num_elem_array);

% This implementation is fast on GPU and for large array
mean_array = convn(convn(convn(input_array, ker1, 'same'), ker2, 'same'), ker3, 'same');
mean_array = mean_array ./ num_elem_array;
% Add abs to avoid numerical error when standard deviation is very small
std_array = sqrt(abs(convn(convn(convn(input_array.^2, ker1, 'same'), ker2, 'same'),...
    ker3, 'same')./ num_elem_array - mean_array.^2));
% Shift the mean back. Ignore the saturation, since the relative correction is tiny.  
mean_array = mean_array - 1;
% Remove the values outside the mask 
% inf_mask = imdilate(compute_mask, strel('cube',3));
switch out_mask_value
    case 'nan'
        std_array(~mask) = nan;
        mean_array(~mask) = nan;
    case 0
        std_array(~mask) = 0;
        mean_array(~mask) = 0;
    case 'inf'
        std_array(~mask) = inf;
        mean_array(~mask) = inf;
end
if isa(input_array, 'gpuArray') && save_to_RAM_Q
    [mean_array, std_array] = gather(mean_array, std_array);
end
    
end
