function bbox = getBbox(origin, sz)
bbox = zeros(8,3);
for i = 1:8
    bbox(i,:) = origin+bin2vec(i,3).*sz;
end
