function [transform, maxC, C_valid] = MaskedTranslationRegistration(fixedImage,movingImage,fixedMask,movingMask,search_range, overlapRatio)

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
% 1. Generalize to 3D masked FFT registration and add searching range limits 
% 2. Accelerate the computation with nearly no lost in accuracy in the
% computation of the correlation array. 
if nargin < 5
    search_range = size(movingImage) - 1; % Search the translation only in +- search_range 
    overlapRatio = 3/10;
elseif nargin < 6
    overlapRatio = 3/10;
end
[C,numberOfOverlapMaskedPixels] = normxcorrn_masked(fixedImage,movingImage,fixedMask,movingMask, search_range);
% fixed_image_size = size(fixedImage);
moving_image_size = size(movingImage);
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
    transform = [(xpeak - search_range(2)) (ypeak - search_range(1))];
elseif ndims(C) == 3
    C_valid = C(valid_min(1):valid_max(1), valid_min(2):valid_max(2), valid_min(3):valid_max(3));    
    [maxC, imax] = max(C_valid(:));
    [ypeak, xpeak, zpeak] = ind2sub(size(C_valid),imax(1));
%     fprintf('Peak position (y, x, z) = (%d, %d, %d)\n', ypeak, xpeak, zpeak);
    transform = [(xpeak - search_range(2)), (ypeak - search_range(1)), (zpeak - search_range(3))];
end
transform = transform';
% The resulting translation can be fed into imtranslate to translate the
% moving image directly. 
% Take the negative of the transform so that it has the correct sign.
% transform = -transform;
end

%% Sub functions
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
for tmp_dim = flip_dimension
    rotatedMovingImage = flip(rotatedMovingImage, tmp_dim);
    rotatedMovingMask = flip(rotatedMovingMask, tmp_dim);
end
clear movingImage movingMask;

% Calculate all of the FFTs that will be needed.
% fixedImageSize = size(fixedImage);
% combinedSize = fixedImageSize + padImageSize - 1;
combinedSize = size(rotatedMovingImage) + padImageSize - 1;
% Find the next largest size that is a multiple of a combination of 2, 3,
% and/or 5.  This makes the FFT calculation much faster.
optimalSize = arrayfun(@FindClosestValidDimension, combinedSize);

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

for ifac = [2 3 5]
    while( rem(n,ifac) == 0 )
        n = n/ifac;
    end
end
end