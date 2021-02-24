function [A11, A12, A13, A22, A23, A33] = fun_gamma_normalized_hessian_matrix3D(A11, gf_sigma, gamma)
% fun_hessian_matrix computes the second order derivative of the
% gaussian-smothed version of data use central difference. 
% Input: 
%   A11: input 3D array, merge variable name with A11 for performance. 
%   gf_size: standard deviation of the gaussian filter. default value is 0
%   gamma: parameter for gamma-normalization
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
if nargin < 2
    gf_sigma = 1;
    gamma = 1;
elseif nargin < 3
    gamma = 1;
end
gamma_normalization_factor = gf_sigma ^ gamma;

if gamma_normalization_factor > 1 + eps
    [A11, A22, A33] = fun_gradient3D(A11.* gamma_normalization_factor,[1,2,3]);
    [A11, A12, A13] = fun_gradient3D(A11.* gamma_normalization_factor,[1,2,3]);
    [A22, A23] = fun_gradient3D(A22.* gamma_normalization_factor,[2,3]);
    A33 = fun_gradient3D(A33.* gamma_normalization_factor,3);
else
    [A11, A22, A33] = fun_gradient3D(A11,[1,2,3]);
    [A11, A12, A13] = fun_gradient3D(A11,[1,2,3]);
    [A22, A23] = fun_gradient3D(A22,[2,3]);
    A33 = fun_gradient3D(A33,3);
end

end