function out_image = fun_stretch_contrast(inputImage, saturate_ptl_low, saturate_ptl_high,image_class, ignore_zeroQ)
% fun_stretch_contrast stretch the contrast of the image array of arbitrary
% dimension by saturating the specified percentile of pixels and linear
% transform
% Input:
%   inputImage: numerical array of arbitrary dimension
%   saturate_ptl_low: float [0, 1], lower percentile for saturation. If not
%   provided, default value is 0.005;
%   saturate_ptl_high: float in [0,1[, higher percentile for saturation. If
%   not provided, the default value is 0.9995;
% Output:
%   contrast-stretched image array
if nargin < 2
    saturate_ptl_low = 0.005;
    saturate_ptl_high = 0.9995;    
    image_class = class(inputImage);
    ignore_zeroQ = false;
elseif nargin < 4
    image_class = class(inputImage);
    ignore_zeroQ = false;
elseif nargin <5
    ignore_zeroQ = false;
end


if (saturate_ptl_low > 1) || (saturate_ptl_high > 1)
    warning('The saturate level should be in [0, 1]. Rescale input saturation range automatically.')
    saturate_ptl_low = saturate_ptl_low / 100;
    saturate_ptl_high = saturate_ptl_high / 100;
end
% For debug
% inputImage = rescale(vascBlock);
% saturate_ptl_low = 0.005;
% saturate_ptl_high = 0.9995;
switch image_class
    case 'uint8'
        conversionFcn = @im2uint8;
    case 'uint16'
        conversionFcn = @im2uint16;
    case 'int16'
        conversionFcn = @im2int16;
    case {'single', 'double', 'float'}
        % Do nothing(?)
end

if ignore_zeroQ
    out_image = inputImage(:);
    out_image(inputImage == 0) = [];
else
    out_image = inputImage(:);
end

if ~isempty(out_image)
    out_image = sort(out_image, 'ascend');

    num_pixels = numel(out_image);
    low_limit = double(out_image(round(max(1, saturate_ptl_low * num_pixels))));
    high_limit = double(out_image(round(max(2, saturate_ptl_high * num_pixels))));
    if low_limit ~= high_limit
        out_image  = max(low_limit, min(high_limit, single(inputImage)));
        out_image = (out_image - low_limit) ./ (high_limit - low_limit);
    else
        out_image = inputImage;
    end
else
    out_image = inputImage;
end


switch image_class
    case {'uint8', 'uint16', 'int16'}
        out_image = conversionFcn(out_image);
end

end