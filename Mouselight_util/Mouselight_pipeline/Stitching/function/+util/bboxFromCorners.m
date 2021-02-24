function [box, ori, sz]= bboxFromCorners(pts)
ori = min(pts);
sz = max(pts)-min(pts);

lims = [ori; ori+sz];
box = zeros(8,3);

for i = 0:7
    inds = round([1+mod(i,2), 1+mod(floor(i/2), 2), 1+mod(floor(i/4), 2)]);
    box(i+1,:) = [lims(inds(1), 1), lims(inds(2), 2), lims(inds(3), 3)];
end
end


