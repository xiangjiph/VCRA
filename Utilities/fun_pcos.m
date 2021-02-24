function pcos = fun_pcos(vec_mat)
% fun_pcos computs pairwise vector projection angle
% vec_mat: N-by-d array, N vectors in d dimensional space 
num_vec = size(vec_mat, 1);

vec_mat = vec_mat';
% Normalization
vec_mat_norm = sum(vec_mat .^2, 1);
vec_mat = vec_mat ./ vec_mat_norm;

pcos = nan(num_vec, num_vec);
for iter_vec = 1 : num_vec
    tmp_vec = vec_mat(:, iter_vec)';
    pcos(:, iter_vec) = tmp_vec * vec_mat;
end
    