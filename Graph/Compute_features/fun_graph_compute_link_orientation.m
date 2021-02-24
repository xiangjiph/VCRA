function vec = fun_graph_compute_link_orientation(sub)
% fun_graph_compute_link_orientation compute the orientation of the link
% connected component by PCA. 
% Input: 
%   sub: N-by-3 numerical array, subscript of the link voxels 
% Output:
%   vec: 3-by-1 array, the normalized eigenvector correspond to the greatest
%   eigenvalue.

if size(sub, 1) < 3
    vec = [nan;nan;nan];
else
    [vec] = pca(sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
end

end