function [X, Y] = fun_vis_bbox_mmxx_to_polygon_vertice2D(bbox_mmxx)

% Each colume: left, front, back, right, bottom, top
polyX_ind = [2, 4, 4, 2];
polyY_ind = [1, 1, 3, 3];

if isvector(bbox_mmxx)
    X = bbox_mmxx(polyX_ind);
    Y = bbox_mmxx(polyY_ind);
else
    % If bbox_mmxx is a N-by-6 array, where N is the number of bounding
    % boxes
   X = bbox_mmxx(:, polyX_ind);
   Y = bbox_mmxx(:, polyY_ind);
   
   X = reshape(X', 4, []);
   Y = reshape(Y', 4, []);
end
end