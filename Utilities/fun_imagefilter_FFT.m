function X = fun_imagefilter_FFT(X,ker, padOption, padDirection)
% fun_imagefilter_FFT compute the convolution between the input image and
% kernel.Before convolution, the image is padded according to the input parameters. 
% The convolution is done in frequency space. 
% This function is compatable with gpuArray

% Input: 
%   X: image, 3 dimensional array. If X is not float, it will be converted to
%   single. 
%   ker: kernel, 3 dimensional array
%   padOption: option for padding the image. Input for padarray
%   padDirection: option for padding the image. Input for padarray
% Output:
%   X: image after convolution
if ~exist('padOption', 'var')
    padOption = 'replicate';
end
if ~exist('padDirection', 'var')
    padDirection = 'both';
end

if ~isfloat(X)
    X = single(X);
end
if ~isfloat(ker)
    ker = single(ker);
end


x_size = size(X);
ker_size = size(ker);
switch padDirection
    case 'both'
        X = padarray(X, ceil((ker_size - 1)/2), 'both', padOption);
%     case {'post', 'pre'}
%         X = padarray(X, ker_size - 1, padDirection, padOption);
end
pad_size = size(X);
% Choose best fft size
optimized_fft_size = find_optimized_size(pad_size);
% Convolution theorem
X = ifftn(fftn(X, optimized_fft_size) .* fftn(ker, optimized_fft_size), 'symmetric');
% Remove the padded boundary 
switch padDirection
    case 'both'
        data_start_idx = 1 + pad_size - x_size;
        data_end_idx = data_start_idx + x_size - 1;
        if length(x_size) == 3
            X = X(data_start_idx(1):data_end_idx(1), data_start_idx(2):data_end_idx(2), data_start_idx(3):data_end_idx(3));
        elseif ismatrix(X)
            X = X(data_start_idx(1):data_end_idx(1), data_start_idx(2):data_end_idx(2));
        end
end
end

function optimized_size = find_optimized_size(array_size)
% find_optimized_size compute the optimized FFT size for the given
% array size. The size in each dimension should be smaller than 10000. 


hams = [8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48 ...
            ,50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128 ...
            ,135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250 ...
            ,256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450 ...
            ,480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729 ...
            ,750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080, 1125 ...
            ,1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536 ...
            ,1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160 ...
            ,2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916 ...
            ,3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840 ...
            ,3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800, 4860, 5000 ...
            ,5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144, 6250, 6400 ...
            ,6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776, 8000, 8100 ...
            ,8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000];
        
% next_fast_len = @(l)(hams(find(hams>l, 1, 'first')));
% optimized_size = arrayfun(@(l)next_fast_len(l),array_size);
ndims = length(array_size);
optimized_size = zeros();
for dim_idx = 1 : ndims
    optimized_size(dim_idx) = hams(find(hams >= array_size(dim_idx), 1, 'first'));
end
end