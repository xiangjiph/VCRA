function output = fun_laplacian_of_gaussain3(sigma, filter_size)
% FUN_LAPLACIAN_OF_GAUSSAIN3 creats a laplacian of gaussian filter in 3
% dimension
% Input: 
%   sigma: floating point scalar, the standard deviation for computing the
%      gaussian filter;
%   filter_size: scalar, the size of the filter. The default value is 6
%       times of the sigma + 1;
% Output: 
%   output: the filter. 

if nargin < 2
    filter_size = 6 * sigma + 1;
end

filter_radius = (filter_size - 1)/2;
[x,y,z] = meshgrid( -filter_radius : filter_radius, -filter_radius : filter_radius, -filter_radius : filter_radius);
exp_part = exp( - (x.^2 + y.^2 + z.^2 ) / (2 * sigma^2) );
exp_part(exp_part < eps * max(exp_part(:))) = 0;
sum_exp_part = sum(exp_part(:));
if sum_exp_part ~= 0
    exp_part = exp_part ./ sum_exp_part;
end

output = - (( x.^2  + y.^2 + z.^ 2 - sigma^2) / sigma^4) .* exp_part;
output = output - sum(output(:))/filter_size^3;

end