function eigenvalues = fun_cc_compute_moment_eigenvalues(cc)
% fun_cc_compute_moment_eigenvalues computes the eigenvalues of the second
% moment matrix of the connected component in three dimension. 
% Input:
%   cc: connected components output by bwconcomp
% Output:
%   e1, e2, e3: double precision vector, correspont to three eigenvalues,
%   orderred as |e1| >= |e2| >= |e3|

[im_200, im_110, im_101, im_020, im_011, im_002] = fun_cc_compute_moments(cc);
[e1, e2, e3] = fun_vectorized_3x3_real_sym_mat_eigenvalues(im_200, im_110, im_101, im_020, im_011, im_002);
eigenvalues = cat(2, e1, e2, e3);
end

