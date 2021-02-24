function varargout = fun_gradient3D(inputArray, dim, method, save_to_ramQ )
% fun_gradient3D compute the gradient along dim. Forward difference is used
% at the boundary, central difference is used for the internal points
% Input: 
%   inputArray: 3D numerical array
%   dim: scalar or vector specifying the dimension(s) along which the
%   gradient is computed. 
%   save_to_ramQ: if the inputArray is a gpuArray and save_to_ramQ is true,
%   gather all the result to RAM.
% Output: 
%   varargout: computed gradient(s) 
if nargin < 3
    save_to_ramQ = false;
    method = 'central';
end

if nargin < 4
    save_to_ramQ = false;
end

ndim = length(dim);
varargout = cell(1,ndim);
input_size = size(inputArray);


for dir_idx = 1 : ndim
    grad = zeros(input_size, 'like', inputArray);
    grad_dir = dim(dir_idx);
    switch grad_dir
        case 1
            switch method
                case 'central'
                    grad(1,:,:) = inputArray(2,:,:) - inputArray(1,:,:);
                    grad(input_size(1),:,:) = inputArray(input_size(1),:,:) - inputArray(input_size(1) - 1, :, :);
                    grad(2:input_size(1)-1, :,:) = (inputArray(3:input_size(1),:,:) - inputArray(1:input_size(1) - 2, :, :))./2;
                case 'intermediate'
                    grad(1:end-1,:,:) = diff(inputArray, 1, 1);
            end
        case 2
            switch method
                case 'central'
                    grad(:,1,:) = inputArray(:,2,:) - inputArray(:,1,:);
                    grad(:,input_size(2),:) = inputArray(:,input_size(2),:) - inputArray(:,input_size(2) - 1, :);
                    grad(:, 2:input_size(2) - 1,:) = (inputArray(:,3:input_size(2),:) - inputArray( :,1:input_size(2) - 2, :))./2;
                case 'intermediate'
                    grad(:, 1:end-1,:) = diff(inputArray, 1, 2);
                    
            end
        case 3
            switch method
                case 'central'
                    grad(:,:,1) = inputArray(:,:,2) - inputArray(:,:,1);
                    grad(:,:,input_size(3)) = inputArray(:,:,input_size(3)) - inputArray(:,:,input_size(3) - 1);
                    grad(:,:,2:input_size(3)-1) = (inputArray(:,:,3:input_size(3)) - inputArray(:,:,1:input_size(3) - 2))./2;
                case 'intermediate'
                    grad(:, :, 1:end-1) = diff(inputArray, 1, 3);
            end
    end
    
    if save_to_ramQ && isa(inputArray, 'gpuArray')
        varargout{dir_idx} = gather(grad);
    else
        varargout{dir_idx} = grad;
    end
    
end

end