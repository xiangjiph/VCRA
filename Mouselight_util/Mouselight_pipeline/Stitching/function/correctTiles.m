function desc = correctTiles(desc,dims)
% flip 
if isempty(desc)
    return;
end
desc(:,1:2) = dims(1:2)+1 - desc(:,1:2);
end