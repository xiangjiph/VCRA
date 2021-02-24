function X =convFFT(X,ker, shape, pad_size, pad_method)
% convFFT do convolution of 1/2/3 dimension based on Fourier Transform
% X can be the image array; Ker is the kernel for convolution. 
% Cast both the bot X and ker into single unless X,ker are double or single. 
% X = I;
% ker = fun_gaussian_kernell([sigma_list(sigma_idx),sigma_list(sigma_idx),sigma_list(sigma_idx)]);
% pad_method = '0';
% shape = 'same';

x_size = size(X);
ker_size = size(ker);
fft_size = x_size + ker_size - 1;

if nargin < 4
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
    next_fast_len = @(l)(hams(find(hams>=l, 1, 'first')));
    pad_size = arrayfun(@(l)next_fast_len(l),fft_size);
    pad_method = '0';
elseif nargin < 5
    pad_method = '0';
end


if ~isfloat(X)
    X = single(X);
end
if ~isfloat(ker)
    ker = single(ker);
end

switch pad_method
    case '0'
        X = ifftn(fftn(X,pad_size) .* fftn(ker,pad_size),'symmetric');
    case 'replicate'
        X = ifftn(fftn(padarray(X,pad_size - x_size, 'post', 'replicate')) .* fftn(ker,pad_size),'symmetric');
    case 'symmetric'
        X = ifftn(fftn(padarray(X,pad_size - x_size, 'post', 'symmetric')) .* fftn(ker,pad_size),'symmetric');
end

switch shape
    case 'same'
        output_size = x_size;
    case 'full'
        output_size = fft_size;
    case 'valid'
        output_size = x_size - ker_size + 1;
end

start_ind = ceil((fft_size - output_size) /2) + 1;
end_ind = start_ind + output_size - 1;
if length(x_size) == 3
    X = X(start_ind(1) : end_ind(1), start_ind(2) : end_ind(2), start_ind(3) : end_ind(3));
elseif isvector(ker)
    % It's rediculous that MATLAB views a vector as a matrix......
    X = X(max(start_ind(:)) : min(end_ind(:)));
elseif ismatrix(ker)
    X = X(start_ind(1) : end_ind(1), start_ind(2): end_ind(2));
end

end