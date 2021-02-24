function output_array = moving_std(input_array, window_size,use_fftQ,minus1Q)
% MOVING_STD computes the standard deviation of the given input_array
% within the given window_size. 
% Input: 
%   input_array: n-dimensional array
%   window_size: size of the window for computing the local standard
%       deviation
%   use_fftQ: if true, use FFT for convolution
%   minus1Q: if true, compute the unbaised estimation of the standard deviation
if nargin < 3
    minus1Q = true;
    use_fftQ = true;
elseif nargin < 4
    minus1Q = true;
end
% Convert the array into single precesion to improve the computation on GPU.
% For the purpose in this project, we don't need double precision. 
if isa(input_array, 'double')
    input_array = single(input_array);
end
% If input window_size is a scalar, assume isotropic window size. 
if (ismatrix(input_array)) && (length(window_size)==1)
    window_size = [window_size, window_size];
elseif (ndims(input_array)>2) && (length(window_size)==1)
    window_size = ones([1,ndims(input_array)]) .* window_size;
end

window_array = ones(window_size,'single');
num_elements = length(window_array(:));

if use_fftQ
    output_array = convFFT_v2(input_array.^2, window_array, 'same') - ...
        ((convFFT_v2(input_array, window_array,'same').^2)./num_elements);
else
    output_array = convn(input_array.^2, window_array, 'same') - ...
        ((convn(input_array, window_array,'same').^2)./num_elements);
end
if minus1Q
    output_array = sqrt(output_array./(num_elements-1));
else
    output_array = sqrt(output_array./(num_elements));
end

end
