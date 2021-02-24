function [mean_array, std_array] = fun_compute_local_meanNstd(std_array, window_size,minus1Q,saveOnMemQ)
% COMPUTE_LOCAL_MEANNSTD computes the mean and standard deviation of the given input_array
% within a mvoing window of size window_size. 
% Input: 
%   std_array: input_array, n-dimensional array
%   window_size: size of the window for computing the local standard
%       deviation
%   minus1Q: if true, compute the unbaised estimation of the standard deviation
%   saveOnMemQ: return the GPU array to memory.
% Output: 
%   mean_array: numerical array, same size as input_array. Mean value of the voxel value in the mvoing window
%   input_array: numerical array. Standard deviation of the voxel value in
%   the moving windows. 

if nargin < 3
    minus1Q = true;
    saveOnMemQ = true;
elseif nargin < 4
    saveOnMemQ = true;
end
% Convert the array into single precesion to improve the computation on GPU.
% Actually we might need double precision some time... When the standard
% deviation is very small and thus the numerical error might be large. 
if ~isa(std_array, 'float')
    std_array = single(std_array);
end
% If input window_size is a scalar, assume isotropic window size. 
if (ismatrix(std_array)) && (length(window_size)==1)
    window_size = [window_size, window_size];
elseif (ndims(std_array)>2) && (length(window_size)==1) %#ok<ISMAT>
    window_size = ones([1,ndims(std_array)]) .* window_size;
end

num_elements = prod(window_size);
mean_array = fun_decomposed_convolution_spatial(std_array, window_size, 'replicate', 'both');
std_array = fun_decomposed_convolution_spatial(std_array.^2, window_size, 'replicate', 'both') - ...
        ((mean_array.^2)./num_elements);
% Correction to the negative numerical error when the variance is small. 
std_array = abs(std_array);
if minus1Q
    std_array = sqrt(std_array./(num_elements-1));
else
    std_array = sqrt(std_array./(num_elements));
end
mean_array = mean_array ./ num_elements;

if isa(mean_array, 'gpuArray') && saveOnMemQ
    mean_array = gather(mean_array);
    std_array = gather(std_array);
end
end



%% Subfunction
function input_array = fun_decomposed_convolution_spatial(input_array,window_size, padValue, padDirection)
x_size = size(input_array);
if nargin < 3
    padValue = 'replicate';
    padDirection = 'both';
elseif nargin < 4
    padDirection = 'both';
end

if ~isfloat(input_array)
    input_array = single(input_array);
end
% Construct convolution kernels in three directions
ker1 = ones(window_size(1),1);
ker2 = ones(1, window_size(2));
ker3 = ones(1, 1, window_size(3));

pad_width = ceil((window_size - 1)/2);
switch padDirection
    case 'both'
        input_array = padarray(input_array, pad_width, 'both', padValue);
    otherwise
        error('Unrecognized pad direction');
end
% Decompose the convolution into three one-dimensional convolution. 
% X = convn(X, window_size, 'same');
input_array = convn(convn(convn(input_array, ker1, 'same'), ker2, 'same'), ker3, 'same');

switch padDirection
    case 'both'
        data_start_idx = 1 + pad_width;
        data_end_idx = data_start_idx + x_size - 1;
        if length(x_size) == 3
            input_array = input_array(data_start_idx(1):data_end_idx(1), data_start_idx(2):data_end_idx(2), data_start_idx(3):data_end_idx(3));
        elseif ismatrix(input_array)
            input_array = input_array(data_start_idx(1):data_end_idx(1), data_start_idx(2):data_end_idx(2));
        end
end
end