function sub = fun_ind2sub(block_size, ind)

if numel(block_size) == 2
    sub = zeros(numel(ind), 2);
    [sub(:,1), sub(:,2)] = ind2sub(block_size, ind);
elseif numel(block_size) == 3
    sub = zeros(numel(ind), 3);
    [sub(:,1), sub(:,2), sub(:,3)] = ind2sub(block_size, ind);
else
    error('Unsupported block size');
end
end