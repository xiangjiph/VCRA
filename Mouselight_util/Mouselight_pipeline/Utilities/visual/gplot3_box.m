function h = gplot3_box(A, xyz, bbox, varargin)
%GPLOT Plot graph (nodes and edges).
%   GPLOT(A, xyz) plots the graph specified by the adjacency matrix,
%   A, and the n-by-3 coordinate array, xyz.
%   
%   GPLOT(A, xyz, linespec) uses line type and color specified in the
%   string LineSpec. See PLOT for possibilities.
%
%   h = GPLOT(A, xyz) returns the a handle to the graph.
%   
%   h = GPLOT(A, xyz, 'LineWidth', 5, ...) also takes arbitrary arguments
%   for line properties
    % If no arguments given, then run buckminster sphere example
    if nargin == 0
        [A, xyz] = bucky;
    end

    % If only one argument given, throw error.
    if nargin == 1
        error('Please provide an adjacency matrix and coordinate array');
    end
    
    if nargin <3
        bbox=[];
    end

    % Returns i and j, lists of connected nodes
    [i,j] = find(A);

    % Extact 
    X = [ xyz(i,1) xyz(j,1)]';
    Y = [ xyz(i,2) xyz(j,2)]';
    Z = [ xyz(i,3) xyz(j,3)]';
    
    % Add NaN values to break between line segments
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
    Z = [Z; NaN(size(i))'];

    % Filter
    if ~isempty(bbox)
        %bbox = [20066+[-100 100] 10298+[-100 100] 2959+[-100 100]];
        it = X>bbox(1)&X<bbox(2)&Y>bbox(3)&Y<bbox(4)&Z>bbox(5)&Z<bbox(6);
        it = all(it(1:2,:));
        X = X(:,it);
        Y = Y(:,it);
        Z = Z(:,it);
    end
    
    % Serialize the x and y data
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    % If only two arguments, then plot as is
    if nargin == 0 || nargin <4
        h = plot3(X, Y, Z);
    end
    
    % If linespec given, then use it
    if nargin >= 4
        if mod(nargin, 2) == 1
            h = plot3(X, Y, Z, varargin{1});
            start = 2;
        else
            h = plot3(X, Y, Z);
            start = 1;
        end
        
        % Now apply the rest of the var string
        if ~isempty(varargin)
            for i=start:2:length(varargin)
                set(h, varargin{i}, varargin{i+1});
            end
        end
        
    end
    
end