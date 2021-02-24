function [outputArgs] = plotRects(rects,fig,view_)
%plotRects plot rectangle array. Helper function for neuron crawler
%
% [OUTPUTARGS] = plotRects(INPUTARGS) Explain usage here
%
% Inputs:
%   rects: [8x3xN] : locations of each corner of BB in 3D
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/09/16 13:39:41 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if nargin<3
    view_ = [20,34];
end
edges = [[1 2];[1 3];[3 4];[4 2];[1 5];[2 6];[5 6];[5 7];[6 8];[3 7];[4 8];[7 8]];
X = [];
C = [];
Col = nan(size(rects,3)*3,3);
if size(rects,3)==1
    Col(:) = 0;
else
    Col(1:3:end,:) = jet(size(rects,3));
    Col(2:3:end,:) = jet(size(rects,3));
    Col(3:3:end,:) = jet(size(rects,3));
end
for i=1:size(edges,1)
    Xtmp =nan(3,size(rects,3)*3);
    tmp1 = squeeze(rects(edges(i,1),:,:));
    tmp2 = squeeze(rects(edges(i,2),:,:));
    Xtmp(:,1:3:end) = tmp1;
    Xtmp(:,2:3:end) = tmp2;
    X = [X Xtmp];
    C = [C;Col];
end
X = X';
if nargin>1
    if fig
        h = figure(fig);
    else
        
    end
else
    h = figure;
end
% set(h,'Position',[-3758 1 1879 2053])
view(view_)
if 0
    plot3(X(:,1),X(:,2),X(:,3)),
else
    for i=1:size(X,1)/3
%         line(X((i-1)*3+1:i*3,1),X((i-1)*3+1:i*3,2),X((i-1)*3+1:i*3,3),'Color',[0 0 0]),
        line(X((i-1)*3+1:i*3,1),X((i-1)*3+1:i*3,2),X((i-1)*3+1:i*3,3),'Color',C(i*3,:)),
    end
end
axis equal tight,
end
