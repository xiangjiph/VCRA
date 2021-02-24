function input_data = fun_downsample_by_block_operation(input_data, block_fun, downsample_rate, padarray_Q)
% fun_downsample_by_block_operation downsamples the input_data by dividing the
% data into blocks of size specified by the downsample_rate, computing the
% maximum value of each blocks to get the downsampled version of the image.
% Input: 
%   input_data: 3D numerical array
%   block_fun: function handle
%   downsample_rate: scalar or 3-by-1 numerical array, should be larger
%   than 1
%   padarray_Q: pad array before downsampling
% Written by Xiang JI on 11/27/2018

if nargin < 4    
    padarray_Q = true;
end

% input_data = image_overview;
% downsample_rate = [2,2,2];
data_size = size(input_data);
data_dim = ndims(input_data);
if data_dim == 2
    data_size(3) = 1;
end

target_data_size = round(data_size ./ downsample_rate);
if isscalar(downsample_rate)
    downsample_rate = ones(1, data_dim) .* downsample_rate;
end

if padarray_Q
    output_data_size = ceil(data_size./downsample_rate);
    pad_array_size = output_data_size .* downsample_rate;
    pad_size = pad_array_size - data_size;
    if any(pad_size)
%         disp('Pad array before downsampling');
        input_data = padarray(input_data, pad_size, 'symmetric', 'post');
        data_size = pad_array_size;
    end
elseif any(mod(data_size, downsample_rate))
    error('Ratio between the input data size and the downsample rate should be integer');
else
    output_data_size = target_data_size;
end

if downsample_rate(1) ~= 1
    input_data = reshape(input_data, downsample_rate(1), prod(data_size)/downsample_rate(1));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(1), data_size(2), data_size(3));
end
if downsample_rate(2) ~= 1
    input_data = permute(input_data, [2,1,3]);
    input_data = reshape(input_data, downsample_rate(2), prod(data_size)/prod(downsample_rate(1:2)));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(2), output_data_size(1), data_size(3));
    input_data = permute(input_data, [2,1,3]);
end
if downsample_rate(3) ~= 1
    input_data = permute(input_data, [3,1,2]);
    input_data = reshape(input_data, downsample_rate(3), prod(output_data_size));
    input_data = block_fun(input_data);
    input_data = reshape(input_data, output_data_size(3), output_data_size(1), output_data_size(2));
    input_data = permute(input_data, [2,3,1]);
end

end