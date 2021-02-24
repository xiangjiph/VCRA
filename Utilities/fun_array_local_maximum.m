function is_local_max_Q = fun_array_local_maximum(A, window_size)
% fun_array_local_maximum computes the local maximum logical array of the
% input array A. 
% Input: 
%   A: 3D numerical array
%   window_size: size of the moving maximum window. 
%
if window_size == 1
    is_local_max_Q = true(size(A));
    return;
end
dt_mov_max = movmax(A, window_size, 1, 'omitnan');
if ~isvector(A)
    num_dim = ndims(A);
    for iter_dim = 2 : num_dim
    dt_mov_max = movmax(dt_mov_max, window_size, iter_dim, 'omitnan');
    end
end
is_local_max_Q = (A >= dt_mov_max);
end