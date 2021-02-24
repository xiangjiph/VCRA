function [maxarg] = findInflection(X)

% make sure data is sorted
if any(X-sort(X,'descend'))
    error('data is not sorted')
end
% first scan
x1 = [1;X(1)];
x2 = [length(X);X(end)];
x0 = [[1:length(X)];[X']];
[~,maxarg] = dist2line(x1,x2,x0);
end
