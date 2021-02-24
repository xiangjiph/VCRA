function vis_implay_in_pair(array1, array2, join_dim)
size1 = size(array1);
size2 = size(array2);
dim1 = ndims(array1);
dim2 = ndims(array2);
if class(array1) ~= class(array2)
    array2 = cast(array2, class(array1));
end

if (dim1 ~= dim2)
    error("The array dimension does not match");
end
max_size = max(size1, size2);
max_size(join_dim) = size1(join_dim) + size2(join_dim);
join_array = ones(max_size, class(array1));
switch join_dim
    case 1
        join_array(1:size1(1),:,:) = array1;
        join_array(size1(1)+1:end,:,:) = array2;
    case 2
        join_array(:, 1:size1(2),:) = array1;
        join_array(:, size1(2)+1:end,:) = array2;
    case 3
        join_array(:,:,1:size1(3)) = array1;
        join_array(:,:,size1(3)+1:end) = array2;
    otherwise
        error("The join dimension should be 1/2/3");
end
implay(rescale(join_array));
end