function [X, Y, Z] = fun_vis_bbox_mmxx_to_polygon_vertice(bbox_mmxx)

% Each colume: left, front, back, right, bottom, top
polyX_ind = [2, 2, 2, 5, 2, 2; ...
    2, 5, 5, 5, 5, 5; ...
    2, 5, 5, 5, 5, 5; ...
    2, 2, 2, 5, 2, 2];
polyY_ind = [1, 4, 1, 1, 1, 1; ...
    1, 4, 1, 1, 1, 1; ...
    4, 4, 1, 4, 4, 4; ...
    4, 4, 1, 4, 4, 4];
polyZ_ind = [3, 3, 3, 3, 6, 3; ...
    6, 3, 3, 6, 6, 3; ...
    6, 6, 6, 6, 6, 3; ...
    3, 6, 6, 3, 6, 3];

if isvector(bbox_mmxx)
    X = bbox_mmxx(polyX_ind);
    Y = bbox_mmxx(polyY_ind);
    Z = bbox_mmxx(polyZ_ind);
else
    % If bbox_mmxx is a N-by-6 array, where N is the number of bounding
    % boxes
   X = bbox_mmxx(:, polyX_ind);
   Y = bbox_mmxx(:, polyY_ind);
   Z = bbox_mmxx(:, polyZ_ind);
   
   X = reshape(X', 4, []);
   Y = reshape(Y', 4, []);
   Z = reshape(Z', 4, []);    
end
end