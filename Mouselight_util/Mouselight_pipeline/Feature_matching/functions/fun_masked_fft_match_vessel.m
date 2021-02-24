function [pixshift_xyz, max_xcorr] = fun_masked_fft_match_vessel(descriptor_1, descriptor_2, pixel_shift_0, iadj, vis_Q)
% fun_masked_fft_match_vessel use masked-fft based cross-correlation to
% determine the translational registration parameter between tile 1 and
% tile 2. 
% Input: 
%   descriptor_1 & descriptor_2: output of vesselDescriptor.m. MATLAB
%   structure including fields: 
%       record.fp_image: filepath to the raw image
%       record.valid_bbox_mmxx: valid part of the image [sub_1_min,
%       sub_2_min, sub_3_min, sub_1_max, sub_2_max, sub_3_max]
%   pixel_shift_0: initial estimation of the pixel shift between two tiles.
%   [x_shift, y_shift, z_shift].
%   iadj: equals 1 or 2 or 3. 1 for x, 2 for y, 3 for z
%   vis_Q: visualize the registration result. 
% Output: 
%   pixshift_xyz: pixel shift between tile 1 and tile 2. [x, y, z]
%   max_xcorr: cross correlation value between two images
% 
% Author: Xiang Ji, UC San Diego, xiangji.ucsd@gmail.com
% Dec 20, 2018
if nargin < 4
    [~, iadj] = max(abs(pixel_shift_0));
end
if nargin < 5
    vis_Q = false;
end

% Parameters - Need to change this number according to microscope? 
empty_pixel_size_xyz = [40, 30, 0]; % The size of the empty region in the raw data from mouselight
switch iadj
    case 1
        search_range_yxz = [10,15,5];
    case 2
        search_range_yxz = [15,10,5];
    case 3
        search_range_yxz = [140, 140, 100];
end
mask_seach_expansion = 0;
% Load images
[pixshift_xyz, max_xcorr] = deal([]);
if isfield(descriptor_1.record, 'fp_image') && isfile(descriptor_1.record.fp_image)
%     fprintf('Loading fixed image stack %s\n', descriptor_1.record.fp_image);
    tile_image_fixed = deployedtiffread(descriptor_1.record.fp_image);
%     fprintf('Fixed image size is (%d, %d, %d)\n', size(tile_image_fixed));
else
    warning('The input descriptor structure does not contains file path to the raw image');
    return;
end
if isfield(descriptor_2.record, 'fp_image') && isfile(descriptor_2.record.fp_image)
%     fprintf('Loading moving image stack %s\n', descriptor_2.record.fp_image);
    tile_image_moving = deployedtiffread(descriptor_2.record.fp_image);
%     fprintf('Moving image size is (%d, %d, %d)\n', size(tile_image_moving));
else
    warning('The input descriptor structure does not contains file path to the raw image');
    return;
end
% Flip the image
tile_image_fixed = flip(flip(tile_image_fixed, 1), 2);
tile_image_moving = flip(flip(tile_image_moving, 1), 2);
% The pixshift is [x_shift, y_shift, z_shift], while the bounding
% box etc are in [y, x, z]. Flip the pixshift for intensity
% registration here.
% if isfield(descriptor_1.record, 'valid_bbox_mmxx')
%     descriptor_valid_bbox_mmxx = descriptor_1.record.valid_bbox_mmxx;
% else
%     descriptor_valid_bbox_mmxx = [8, 45, 9, 1529, 986, 251]; % For 2018-08-15 dataset
% %     descriptor_valid_bbox_mmxx = [90, 120, 40, 1529, 900, 251];
% end

pixshift_yxz0 = pixel_shift_0([2,1,3]);
empty_pixel_shift = [0,0,0];
empty_pixel_shift(iadj) = empty_pixel_size_xyz(iadj);
empty_pixel_shift_yxz = empty_pixel_shift([2,1,3]);
tile_size_yxz = size(tile_image_fixed);

overlap_bbox_1_mmxx = [1, 1, 1, descriptor_1.record.raw_image_size];
overlap_bbox_2_mmxx = [1, 1, 1, descriptor_2.record.raw_image_size];
if pixel_shift_0(iadj) > 0
    overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), pixshift_yxz0 + empty_pixel_shift_yxz) - mask_seach_expansion;
    overlap_bbox_2_mmxx(4:6) = tile_size_yxz - (pixshift_yxz0) + mask_seach_expansion;
elseif pixel_shift_0(iadj) < 0
    overlap_bbox_2_mmxx(1:3) = max(overlap_bbox_2_mmxx(1:3),  - (pixshift_yxz0 - empty_pixel_shift_yxz)) - mask_seach_expansion;
    overlap_bbox_1_mmxx(4:6) = tile_size_yxz + ( pixshift_yxz0 ) + mask_seach_expansion;
end
overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), descriptor_1.record.valid_bbox_mmxx(1:3));
overlap_bbox_2_mmxx(1:3) = max(overlap_bbox_2_mmxx(1:3), descriptor_2.record.valid_bbox_mmxx(1:3));
overlap_bbox_1_mmxx(4:6) = min(overlap_bbox_1_mmxx(4:6), descriptor_1.record.valid_bbox_mmxx(4:6));
overlap_bbox_2_mmxx(4:6) = min(overlap_bbox_2_mmxx(4:6), descriptor_2.record.valid_bbox_mmxx(4:6));

if any(overlap_bbox_1_mmxx(1:3) > overlap_bbox_1_mmxx(4:6))
%     warning('Bounding box for the overlapping region 1 has negative size');
    return;
end
if any(overlap_bbox_2_mmxx(1:3) > overlap_bbox_2_mmxx(4:6))
%     warning('Bounding box for the overlapping region 1 has negative size');
    return;
end
overlap_bbox_1_mmxx(1:3) = min(overlap_bbox_1_mmxx(1:3), overlap_bbox_1_mmxx(4:6));
overlap_bbox_2_mmxx(1:3) = min(overlap_bbox_2_mmxx(1:3), overlap_bbox_2_mmxx(4:6));

overlap_bbox_1_mmll = overlap_bbox_1_mmxx;
overlap_bbox_1_mmll(4:6) = overlap_bbox_1_mmxx(4:6) - overlap_bbox_1_mmxx(1:3) + 1;
overlap_bbox_2_mmll = overlap_bbox_2_mmxx;
overlap_bbox_2_mmll(4:6) = overlap_bbox_2_mmxx(4:6) - overlap_bbox_2_mmxx(1:3) + 1;
% 3D Masked FFT
test_image_fixed = crop_bbox3(tile_image_fixed, overlap_bbox_1_mmll, 'default');
test_image_moving = crop_bbox3(tile_image_moving, overlap_bbox_2_mmll, 'default');
% Threshold for generating the mask 
est_int_th = 1.5e4;
% disp('3D Masked FFT');
% [translation_xyz, max_xcorr, ~] = MaskedTranslationRegistration(test_image_fixed, test_image_moving, ...
%     test_image_fixed > est_int_th , test_image_moving > est_int_th, search_range_yxz);
% Downsample to reduce computation load
test_image_fixed_ds = imresize3(test_image_fixed, 0.5);
test_image_moving_ds = imresize3(test_image_moving, 0.5);
[translation_xyz_ds, max_xcorr, ~] = MaskedTranslationRegistration(test_image_fixed_ds, test_image_moving_ds, ...
    test_image_fixed_ds > est_int_th , test_image_moving_ds > est_int_th, ceil(search_range_yxz/2));
translation_xyz = round(translation_xyz_ds * 2);
pixshift_xyz = overlap_bbox_1_mmll([2,1,3]) - overlap_bbox_2_mmll([2,1,3]) + translation_xyz';
if vis_Q
    vis_pixshift_xyz = pixshift_xyz;
    vis_translation = vis_pixshift_xyz - overlap_bbox_1_mmll([2,1,3]) + overlap_bbox_2_mmll([2,1,3]);
    [test_image_2_moved]= imtranslate(test_image_moving, vis_translation);
%     vis_sec = 10;
%     vis_image_1 = test_image_fixed(:, :, vis_sec);
%     vis_image_2 = test_image_moving(:, :, vis_sec - vis_translation(3));
%     vis_image_2_moved = test_image_2_moved(:, :, vis_sec);
%     Max projection 
    vis_image_1 = max(test_image_fixed, [], 3);
    vis_image_2 = max(test_image_moving, [], 3);
    vis_image_2_moved = max(test_image_2_moved, [], 3);
    figure;
    subplot(1,4,1);
    imshow(vis_image_1);
    title('Tile 1 max projection');
    %         title('Section from tile 1');
    subplot(1,4,2)
    imshow(vis_image_2);
    %         title('Section from tile 2');
    title('Tile 2 max projection');
    subplot(1,4,3)
    imshow(vis_image_2_moved);
    title('Translated tile 2 max projection');
    %         title('Translated section from tile 2');
    subplot(1,4,4)
    imshowpair(vis_image_1, vis_image_2_moved);
    title(sprintf('Image overlap: pixel shift (%d, %d, %d)', vis_pixshift_xyz));
end
end
%% Sub function
function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
% 
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end
warning on
end
%% crop bounding box
function output = crop_bbox3(data, bbox_parameters, bbox_order)
% CROP_BBOX3 crops part of the array DATA according to the given bounding
% box parameters. 
% default bbox_parameters = [ul1, ul2, ul3, l1, l2, l3]
% matlab's regionpros3 output bbox = [ul2, ul1, ul3, l2, l1, l3]
if nargin < 3
    bbox_order = 'default';
    warning('bbox_parameters order not specify. Option: default/ regionprop');
end
if ~iscell(bbox_order)
    bbox_parameters = num2cell(round(bbox_parameters));
end
switch bbox_order
    case {'default'}
        if length(bbox_parameters) <=4
            [ul1, ul2, ul3, l1] = bbox_parameters{:};
            l2 = l1;
            l3 = l1;
        else
            [ul1, ul2, ul3, l1, l2, l3] = bbox_parameters{:};
        end
    case {'regionprop'}
        if length(bbox_parameters) <=4
            [ul2, ul1, ul3, l1] = bbox_parameters{:};
            l2 = l1;
            l3 = l1;
        else
            [ul2, ul1, ul3, l2, l1, l3] = bbox_parameters{:};
        end
end
output = data(ul1:ul1+l1-1, ul2:ul2+l2-1, ul3:ul3+l3-1);
end
%% Masked fft 
function [transform_xyz, maxC, C_valid] = MaskedTranslationRegistration(fixedImage,movingImage,fixedMask,movingMask,search_range, overlapRatio)

% [transform,maxC,C,numberOfOverlapMaskedPixels] =
% MaskedTranslationRegistration(fixedImage,movingImage,fixedMask,movingMask,overlapRatio) 
%   Masked FFT normalized cross-correlation registration of movingImage and
%   fixedImage under masks movingMask and fixedMask.
%   movingMask and fixedMask should consist of only 1s and 0s, where 1
%   indicates locations of useful information in the corresponding image,
%   and 0 indicates locations that should be masked (ignored).
%   fixedImage and movingImage need not be the same size, but fixedMask
%   must be the same size as fixedImage, and movingMask must be the same
%   size as movingImage.
%   If a mask is not needed for either the fixedImage or the movingImage,
%   the fixedMask and/or movingMask can be set to an image of all ones of
%   the same size as the corresponding fixedImage and/or movingImage.
%   The optional overlapRatio specifies the number of pixels needed in the
%   overlap region for meaningful results.  It is specified as a ratio of the
%   maximum number of overlap pixels.  Regions in the resulting correlation
%   image that have fewer than this number of pixels will be set to 0.   
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com

% Modified by Xiang Ji, UC San Diego, xiangji.ucsd@gmail.com
% 1. Generalize to 3D masked FFT registration and add restriction by
% allowing choosing array padding size
% 2. Accelerate the computation with nearly no lost in accuracy in the
% computation of the correlation array. 
if nargin < 5
    search_range = size(movingImage) - 1; % Search the translation only in +- search_range 
    overlapRatio = 3/10;
elseif nargin < 6
    overlapRatio = 3/10;
end
moving_image_size = size(movingImage);
if numel(moving_image_size) == 2
    moving_image_size = [moving_image_size, 1];
end
search_range = min(search_range, moving_image_size - 1);
% disp('Normalized corss correlation');
[C,numberOfOverlapMaskedPixels] = normxcorrn_masked(fixedImage,movingImage,fixedMask,movingMask, search_range);
% fixed_image_size = size(fixedImage);
valid_min = moving_image_size - search_range;
valid_max = min(size(C), moving_image_size + search_range);
% imageSize = size(movingImage);

% Mask the borders;
numberOfPixelsThreshold = overlapRatio * max(numberOfOverlapMaskedPixels(:));
C(numberOfOverlapMaskedPixels < numberOfPixelsThreshold) = 0;
% Take the valid part of C:
if ismatrix(C)
    C_valid = C(valid_min(1):end, valid_min(2):end);    
    [maxC, imax] = max(C_valid(:));
    [ypeak, xpeak] = ind2sub(size(C_valid),imax(1));
    transform_xyz = [(xpeak - search_range(2)) (ypeak - search_range(1))];
elseif ndims(C) == 3
    C_valid = C(valid_min(1):valid_max(1), valid_min(2):valid_max(2), valid_min(3):valid_max(3));    
    [maxC, imax] = max(C_valid(:));
    [ypeak, xpeak, zpeak] = ind2sub(size(C_valid),imax(1));
%     fprintf('Peak position (y, x, z) = (%d, %d, %d)\n', ypeak, xpeak, zpeak);
    transform_xyz = [(xpeak - search_range(2)), (ypeak - search_range(1)), (zpeak - search_range(3))];
end
transform_xyz = transform_xyz';
% The resulting translation can be fed into imtranslate to translate the
% moving image directly. 
% Take the negative of the transform so that it has the correct sign.
% transform = -transform;
end
%% Masked fft
function [C,numberOfOverlapMaskedPixels] = normxcorrn_masked(fixedImage, movingImage, fixedMask, movingMask, padImageSize)

% [C,numberOfOverlapMaskedPixels] =
% normxcorr2_masked(fixedImage, movingImage, fixedMask, movingMask)
%   Masked normalized two-dimensional cross-correlation.
%   This function calculates the Masked NCC using FFTs
%   instead of spatial correlation.  It is therefore much faster for
%   larger structuring elements.
%   The masked normalized cross-correlation of fixedImage and
%   movingImage are calculated using masks fixedMask and movingMask.
%   fixedMask and movingMask should consist of only 1s and 0s, where 1
%   indicates locations of useful information in the corresponding image,
%   and 0 indicates locations that should be masked (ignored).
%   fixedImage and movingImage need not be the same size, but fixedMask
%   must be the same size as fixedImage, and movingMask must be the same
%   size as movingImage.
%   The resulting matrix C contains correlation coefficients and its values
%   range from -1.0 to 1.0. 
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
if nargin < 5
    padImageSize = size(movingImage) - 1;
end
assert(all(padImageSize > 0 ), 'Negative pad image size');
fixedImage = shiftData(fixedImage);
movingImage = shiftData(movingImage);
if islogical(fixedMask) && islogical(movingMask)
    fixedMask = single(fixedMask);
    movingMask = single(movingMask);
else
    % Ensure that the masks consist of only 0s and 1s.  Anything less than or
    % equal to 0 gets set to 0, and everything else gets set to 1.
    fixedMask(fixedMask <= 0 ) = 0;
    fixedMask(fixedMask > 0) = 1;
    movingMask(movingMask <= 0) = 0;
    movingMask(movingMask > 0) = 1;
end
disp('Mask the images');
% The fixed and moving images need to be masked for the equations below to
% work correctly.
fixedImage = fixedImage .* fixedMask;
movingImage = movingImage .* movingMask;

% Flip the moving image and mask in both dimensions so that its correlation
% can be more easily handled. - Convert the correlation to the convolution
% of the flipped image array 
flip_dimension = find(size(movingImage)~=1);
rotatedMovingImage = movingImage;
rotatedMovingMask = movingMask;
% disp('Flip dimensions');
for tmp_dim = flip_dimension
    rotatedMovingImage = flip(rotatedMovingImage, tmp_dim);
    rotatedMovingMask = flip(rotatedMovingMask, tmp_dim);
end
clear movingImage movingMask;

% Calculate all of the FFTs that will be needed.
% fixedImageSize = size(fixedImage);
% combinedSize = fixedImageSize + padImageSize - 1;
rotatedMovingImage_size = size(rotatedMovingImage);
if numel(rotatedMovingImage_size) == 2
    rotatedMovingImage_size = [rotatedMovingImage_size, 1];
end
combinedSize = rotatedMovingImage_size + padImageSize - 1;
% Find the next largest size that is a multiple of a combination of 2, 3,
% and/or 5.  This makes the FFT calculation much faster.
assert(all(combinedSize), 'Exist nonpositve dimension');
optimalSize = arrayfun(@FindClosestValidDimension, combinedSize);
% disp('Compute FFTs');
% Only 6 FFTs are needed.
fixedFFT = fftn(fixedImage,optimalSize);
rotatedMovingFFT = fftn(rotatedMovingImage,optimalSize);
fixedMaskFFT = fftn(fixedMask,optimalSize);
rotatedMovingMaskFFT = fftn(rotatedMovingMask,optimalSize);

% Only 6 IFFTs are needed.
% Compute and save these results rather than computing them multiple times.
numberOfOverlapMaskedPixels = real(ifftn(rotatedMovingMaskFFT.*fixedMaskFFT));
numberOfOverlapMaskedPixels = round(numberOfOverlapMaskedPixels);
numberOfOverlapMaskedPixels = max(numberOfOverlapMaskedPixels,eps);
maskCorrelatedFixedFFT = real(ifftn(rotatedMovingMaskFFT.*fixedFFT));
maskCorrelatedRotatedMovingFFT = real(ifftn(fixedMaskFFT.*rotatedMovingFFT));

numerator = real(ifftn(rotatedMovingFFT.*fixedFFT)) - ...
    maskCorrelatedFixedFFT .* maskCorrelatedRotatedMovingFFT ./ numberOfOverlapMaskedPixels;
clear rotatedMovingFFT fixedFFT

fixedSquaredFFT = fftn(fixedImage.*fixedImage,optimalSize);
fixedDenom = real(ifftn(rotatedMovingMaskFFT.*fixedSquaredFFT)) - ...
    maskCorrelatedFixedFFT.^2 ./ numberOfOverlapMaskedPixels;
fixedDenom = max(fixedDenom,0);
clear rotatedMovingMaskFFT fixedSquaredFFT maskCorrelatedFixedFFT;

rotatedMovingSquaredFFT = fftn(rotatedMovingImage.*rotatedMovingImage,optimalSize);
movingDenom = real(ifftn(fixedMaskFFT.*rotatedMovingSquaredFFT)) - ...
    maskCorrelatedRotatedMovingFFT.^2 ./ numberOfOverlapMaskedPixels;
movingDenom = max(movingDenom,0);
clear fixedMaskFFT rotatedMovingSquaredFFT maskCorrelatedRotatedMovingFFT;    

denom = sqrt(fixedDenom .* movingDenom);
clear fixedDenom movingDenom;
disp('Compute the cross correlation');
tol = 1000*eps( max(abs(denom(:))) );
C = min(max((numerator ./ (denom + eps)) .* (denom > tol), -1), 1);
% Crop out the correct size.
if ismatrix(C)
    C = C(1:combinedSize(1),1:combinedSize(2));
    numberOfOverlapMaskedPixels = numberOfOverlapMaskedPixels(1:combinedSize(1),1:combinedSize(2));
elseif ndims(C) == 3
    C = C(1:combinedSize(1),1:combinedSize(2), 1:combinedSize(3));
    numberOfOverlapMaskedPixels = numberOfOverlapMaskedPixels(1:combinedSize(1), 1:combinedSize(2), 1:combinedSize(3));
end

end
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
function B = shiftData(A)

B = double(A);

is_unsigned = isa(A,'uint8') || isa(A,'uint16') || isa(A,'uint32');
if ~is_unsigned
    
    min_B = min(B(:)); 
    
    if min_B < 0
        B = B - min_B;
    end
    
end
end
%-----------------------------------------------------------------------------
function [newNumber] = FindClosestValidDimension(n)

% Find the closest valid dimension above the desired dimension.  This
% will be a combination of 2s, 3s, and 5s.

% Incrementally add 1 to the size until
% we reach a size that can be properly factored.
if n <= 0
    newNumber = 0;
    return;
end
newNumber = n;
result = 0;
newNumber = newNumber - 1;
while( result ~= 1 )
    newNumber = newNumber + 1;
    result = FactorizeNumber(newNumber);
end
end
%-----------------------------------------------------------------------------
function [n] = FactorizeNumber(n)

if n < 0
    error('n should be non negative integer');
elseif n == 0
    return;
end
for ifac = [2 3 5]
    while( rem(n,ifac) == 0 )
        n = n/ifac;
    end
end
end