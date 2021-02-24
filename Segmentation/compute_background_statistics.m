function [mean_array, std_array] = compute_background_statistics(input_array, sampling_array,use_fftQ,minus1Q,saveOnMemQ)
% MOVING_STD computes the standard deviation of the given input_array
% within the given window_size. 
% Input: 
%   input_array: n-dimensional array
%   sampling_array: logical array of the window for computing the local standard
%       deviation
%   use_fftQ: if true, use FFT for convolution
%   minus1Q: if true, compute the unbaised estimation of the standard deviation
% This implement can be run on both CPU and GPU. The final result is always returned to the memory. 
if nargin < 3
    minus1Q = true;
    use_fftQ = true;
    saveOnMemQ = true;
elseif nargin < 4
    minus1Q = true;
    saveOnMemQ = true;
end
% Convert the array into single precesion to improve the computation on GPU.
% For the purpose in this project, we don't need double precision. 
if ~isa(input_array, 'single')
    input_array = single(input_array);
end

if isa(sampling_array, 'logical')
    num_elements = nnz(sampling_array);
else 
    num_elements = sum(sampling_array(:));
end

if use_fftQ
    mean_array = convFFT_v2(input_array, sampling_array,'same');
    std_array = convFFT_v2(input_array.^2, sampling_array, 'same') - ...
        ((mean_array.^2)./num_elements);
else
    mean_array = convn(input_array, sampling_array,'same');
    std_array = convn(input_array.^2, sampling_array, 'same') - ...
        ((mean_array.^2)./num_elements);
end


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
