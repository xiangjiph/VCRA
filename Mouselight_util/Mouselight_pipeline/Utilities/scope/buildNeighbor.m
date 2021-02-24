function [neighbors] = buildNeighbor(grididx)
%BUILDNEIGHBOR Given grid location, builds an connectivity matrix and
%neighbor edge graph
% 
% [NEIGHBORS] = BUILDNEIGHBOR(GRIDIDX) 
% 
% Inputs: [grididx]: Nx3: tile subscripts locations
% 
% Outputs: [neighbors]: Nx7: neighbors list in [idx left top right bottom
% above below]:[id -x -y +x +y -z +z] format
% 
% Examples: 
% [aa,bb,cc] = ndgrid([-1:1],[-1:1],[-1:1]);
% grididx = 5+[aa(:),bb(:),cc(:)];
% grididx(1:3:end,:) = [];
% [neighbors] = buildNeighbor(grididx)
% 
% See also: 

% $Author: base $	$Date: 2016/08/19 11:01:15 $	$Revision: 0.1 $
% Copyright: HHMI 2016

dims = max(grididx);
N = size(grididx,1);
neighbors = nan(N,7); %[idx left top right bottom above below]:[id -x -y +x +y -z +z]
sixneig = [
    [-1 0 0]; % left (-x)
    [0 -1 0]; % top (-y)
    [1 0 0]; % right (+x)
    [0 1 0]; % bottom (+y)
    [0 0 -1]; % above (-z)
    [0 0 1]]; % below (+z)
indsgrid = sub2ind(dims,grididx(:,1),grididx(:,2),grididx(:,3));
indIM = NaN(dims);
indIM(indsgrid) = 1:N;
for idxN = 1:N
    sub = grididx(idxN,:);
    % check 6 neighbor
    set = nan(6,1);
    sub_ = nan(6,3);
    for ix = 1:6
        sub_(ix,:) = sub + sixneig(ix,:);
    end
    validsubs = all(sub_>0&sub_<=ones(6,1)*dims,2);
    sub__ = sub_(validsubs,:);
    set(validsubs) = indIM(sub2ind(dims,sub__(:,1),sub__(:,2),sub__(:,3)));
    neighbors(idxN,:) = [idxN;set];
end
end
