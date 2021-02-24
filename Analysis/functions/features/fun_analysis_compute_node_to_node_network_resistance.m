function resistance = fun_analysis_compute_node_to_node_network_resistance(conductance_mat, gpuQ)
% fun_analysis_compute_node_to_node_network_resistance computes the network
% two-point resistance using the method described in: 
%   F Y Wu 2004 J. Phys. A: Math. Gen Theory of resistor networks: the two
%   point resistance. 
% Input: 
%   conductance_mat: numerical matrix, graph laplacian of the resistor
%   network. L_{i, j} = \sum_j c_{i, j} - c_{i, j}
%   gpuQ: logical scalar. If true, use GPU for computing two point
%   resistance, but not for diagonalization.
% Output: 
%   resistance: numerical matrix of the same size as conductance_mat.
%   r_{ij} is the resistance between node i and node j in the network. This
%   resistance matrix is full and does not represent the direct
%   connectivity. 
% Note: 
% 1. This implementation actually does redundant computation, since the
%   resistance matrix is symmetric. However, by circshift, this
%   implementation is much faster than nested for loop. 
% 2. This function only apply for the conductance matrix that corresponse
% to one graph connected component. The number of zero eigenvalue equals
% the number of graph connected components. 
% Implemented by Xiang Ji on 09/07/2019

assert(issymmetric(conductance_mat));
if nargin < 2
    gpuQ = false;
end

if issparse(conductance_mat)
    conductance_mat = full(conductance_mat);
end

num_vec = size(conductance_mat, 1);
% Diagonalize the conductance matrix. It seems that diagonalizing matrix on
% GPU does not have speed gain. 
tmp_tic = tic;
[eig_vec, eig_val] = eig(conductance_mat);
fprintf('Finish diagonalizing the conductance matrix. Elapse time is %f seconds\n', toc(tmp_tic));
% Take the eigenvalues into a vector
eig_val = eig_val(eye(size(eig_val), 'logical'));
% Assume only the first element is 0 (up to numerical error) and all the
% other eigenvalues are positive. If the conductance matrix have more than
% 1 connected component, this assertion would fail (?). 
assert(issorted(eig_val, 'ascend') && eig_val(2) > 0);
% Remove the first eigen vector and eigen value. 
eig_vec = eig_vec(:, 2:end)';
eig_val_inv = 1 ./ eig_val(2 : end)';
if gpuQ
    eig_vec = gpuArray(eig_vec);
    eig_val_inv = gpuArray(eig_val_inv);
    resistance = nan([num_vec, num_vec], 'gpuArray');
else
    resistance = nan([num_vec, num_vec]);
end
ind_vec = 1 : num_vec;
% For 10000 x 10000 matrix, normally it would take about 3 minutes on
% 1080Ti. 
tmp_tic = tic;
for iter_shift = 1 : num_vec
    ind_vec_shift = circshift(ind_vec, iter_shift);
    tmp_valid_ind = sub2ind([num_vec, num_vec], ind_vec, ind_vec_shift);
    tmp_r_value = eig_val_inv * ((eig_vec - circshift(eig_vec, iter_shift, 2)) .^ 2);
    resistance(tmp_valid_ind) = tmp_r_value;
end
resistance = gather(resistance);
fprintf('Finish computing the network resistance. Elapse time is %f seconds\n', toc(tmp_tic));
assert(issymmetric(resistance));
end
%% For testing
% test_c1 = 1/10;
% test_c2 = 1/30;
% test_c = [2*test_c1, -test_c1, 0, -test_c1;...
%     -test_c1, 2 * test_c1 + test_c2, -test_c1, -test_c2;...
%     0, -test_c1, 2 * test_c1, -test_c1; ...
%     -test_c1, -test_c2, -test_c1, 2 * test_c1 + test_c2];
% test_r = fun_analysis_compute_node_to_node_network_resistance(test_c);
% assert(test_r(1, 3) == 1/test_c1);
my_resist = fun_simulation_blood_flow_get_vessel_resistance(strandRadiusList(:), strandLengthList(:), 'um');