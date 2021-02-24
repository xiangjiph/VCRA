function G = rmedgeBig(G, s, t)
%RMEDGE Remove edges from a large graph. This is a fix for matlab's
%built in rmedge function that is memory efficient for large scale graphs

% % Determine the indices of the edges to be removed.
% if nargin == 2
%   ind = s;
% else
%   ind = findedge(G, s, t);
%   ind(ind == 0) = [];
% end
% % Convert the edge indices to pairs of Node IDs.
% [s, t] = findedge(G, ind);
% Remove edges from the graph.
nn = numnodes(G);
A = G.adjacency;
Amask = sparse([s;t],[t;s],1,nn,nn)>0;
G = graph(A-A.*Amask);
end