function corners = projectCorners(tform, dims)

inpts = zeros(8,3);
corners = zeros(8,3);
for i = 1:8
    inpts(i,:) = bin2vec(i, 3).*dims(1:3);
    temp = [inpts(i,:) 0 1]*tform';
    corners(i,:) = round(temp(1:3));
end



