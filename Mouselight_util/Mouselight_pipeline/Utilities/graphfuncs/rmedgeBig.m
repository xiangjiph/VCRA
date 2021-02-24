function G = rmedgeBig(G, s, t)
%RMEDGE Remove edges from a large graph. This is a fix for matlab's
%built in rmedge function that is memory efficient for large scale graphs
nn = numnodes(G);
A = G.adjacency;
Amask = sparse([s;t],[t;s],1,nn,nn)>0;
G = graph(A-A.*Amask);
end