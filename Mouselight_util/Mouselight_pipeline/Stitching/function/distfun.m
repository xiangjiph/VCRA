function  d = distfun(ZI,ZJ,W)
if nargin<3
    W = [1 2.5 100000];
%     W = [1 2 inf];
end
dvec = ones(size(ZJ,1),1)*ZI-ZJ;
d = sqrt(sum((dvec*diag(W)).*dvec,2));

