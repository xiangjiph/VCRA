function [transform,maxC,C,numberOfOverlapMaskedPixels] = MaskedTranslationRegistration(fixedImage,movingImage,fixedMask,movingMask,overlapRatio)

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
%

if( nargin < 5 )
    overlapRatio = 3/10;
end
% [C,numberOfOverlapMaskedPixels] = normxcorr2_masked(fixedImage,movingImage,fixedMask,movingMask);

[C,numberOfOverlapMaskedPixels] = normxcorrn_masked(fixedImage,movingImage,fixedMask,movingMask);
imageSize = size(movingImage);

% Mask the borders;
numberOfPixelsThreshold = overlapRatio * max(numberOfOverlapMaskedPixels(:));
C(numberOfOverlapMaskedPixels < numberOfPixelsThreshold) = 0;

[maxC, imax] = max(C(:));
if ismatrix(C)
    [ypeak, xpeak] = ind2sub(size(C),imax(1));
    transform = [(xpeak-imageSize(2)) (ypeak-imageSize(1))];
elseif ndims(C) == 3
    [ypeak, xpeak, zpeak] = ind2sub(size(C),imax(1));
    transform = [(xpeak-imageSize(2)), (ypeak-imageSize(1)), (zpeak-imageSize(3))];
end
transform = transform';
% The resulting translation can be fed into imtranslate to translate the
% moving image directly. 
% Take the negative of the transform so that it has the correct sign.
% transform = -transform;
end

%% Sub functions
function [C,numberOfOverlapMaskedPixels] = normxcorrn_masked(varargin)

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
% 
% Modified by Xiang Ji (UC San Diego) on Dec 13
% 1. Generalize to 3D masked FFT registration. It's extremely
% memory expensive. To align two 1532x1024x251 images stack in single
% precision, it takes up to 200+ GB of memory due to the array padding and
% complex array. The memory requirement can be reduced by specifying the
% registration direction or search range( to be implemented )
% 2. Accelerate the computation with nearly no lost in accuracy in the
% computation of the correlation array. 

[fixedImage, movingImage, fixedMask, movingMask] = ParseInputs(varargin{:});
clear varargin;

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
% can be more easily handled.
flip_dimension = find(size(movingImage)~=1);
rotatedMovingImage = movingImage;
rotatedMovingMask = movingMask;
for tmp_dim = flip_dimension
    rotatedMovingImage = flip(rotatedMovingImage, tmp_dim);
    rotatedMovingMask = flip(rotatedMovingMask, tmp_dim);
end
% rotatedMovingImage = rot90(movingImage,2);
% rotatedMovingMask = rot90(movingMask,2);
clear movingImage movingMask;

% Calculate all of the FFTs that will be needed.
fixedImageSize = size(fixedImage);
movingImageSize = size(rotatedMovingImage);
combinedSize = fixedImageSize + movingImageSize - 1;
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
function [fixedImage, movingImage, fixedMask, movingMask] = ParseInputs(varargin)

narginchk(4,4)

fixedImage = varargin{1};
movingImage = varargin{2};
fixedMask = varargin{3};
movingMask = varargin{4};

validateattributes(fixedImage,{'logical','numeric'},{'real','nonsparse','finite'},mfilename,'fixedImage',1)
validateattributes(movingImage,{'logical','numeric'},{'real','nonsparse','finite'},mfilename,'movingImage',2)
validateattributes(fixedMask,{'logical','numeric'},{'real','nonsparse','finite'},mfilename,'fixedMask',3)
validateattributes(movingMask,{'logical','numeric'},{'real','nonsparse','finite'},mfilename,'movingMask',4)

% If either fixedImage or movingImage has a minimum value which is negative, we
% need to shift the array so all values are positive to ensure numerically
% robust results for the normalized cross-correlation.
fixedImage = shiftData(fixedImage);
movingImage = shiftData(movingImage);
end
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