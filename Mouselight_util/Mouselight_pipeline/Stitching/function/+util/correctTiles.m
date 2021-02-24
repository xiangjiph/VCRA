function desc = correctTiles(desc,dims)
% flip 
if nargin<2
    dims = [1024 1536 251];
end
if isempty(desc)
    return;
end
desc(:,1:2) = dims(1:2)+1 - desc(:,1:2);
end