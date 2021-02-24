function [A11, A12, A13, A22, A23, A33] = fun_hessian_matrix3D(A11, gf_size, gatherQ, memClass)
% fun_hessian_matrix computes the second order derivative of the
% gaussian-smothed version of data use central difference. 
% Input: 
%   A11: input 3D array, merge variable name with A11 for performance. 
%   gf_size: standard deviation of the gaussian filter. default value is 0
% Output: 
%   hessian_matrix: structure with fields:
%       gaussian_std: size of the gaussian filter
%       A11, A12,A13, A22, A23, A33: six matrix element of the hessian
%       matrix(which is symmetric), correspond to Ixx, Ixy, Ixz, Iyy, Iyz,
%       Izz
% Note:
%   The data transfer from the GPU to the CPU is the main performance
%   bottleneck for this function. However, computing graident on CPU is
%   much more time comsuming than on GPU. The overall performance is better
%   if the computation is done on GPU and then transfer back to CPU. 
if isa(A11, 'gpuArray')
    data_class = classUnderlying(A11);
else
    data_class = class(A11);
end

if ~ismember(data_class, {'double', 'single', 'float'})
    warning('The inputdata is of class %s. The output data will also be %s', data_class, data_class)
end
if nargin < 2
    gf_size = 0;
    gatherQ = true;
    memClass = 'double';
elseif nargin < 3
    gatherQ = true;
    memClass = 'double';
elseif nargin < 4
    memClass = 'double';
end

if gf_size > 0
    A11 = imgaussfilt3(A11, gf_size);
end

[A11, A22, A33] = fun_gradient3D(A11,[1,2,3]);
[A11, A12, A13] = fun_gradient3D(A11,[1,2,3]);
[A22, A23] = fun_gradient3D(A22,[2,3]);
A33 = fun_gradient3D(A33,3);

if isa(A11, 'gpuArray') && gatherQ
    A11 = cast(gather(A11), memClass);
    A12 = cast(gather(A12), memClass);
    A13 = cast(gather(A13), memClass);
    A22 = cast(gather(A22), memClass);
    A23 = cast(gather(A23), memClass);
    A33 = cast(gather(A33), memClass);
end

end