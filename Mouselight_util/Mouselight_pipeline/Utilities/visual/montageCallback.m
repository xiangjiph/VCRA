function output_txt = montageCallback(obj,event_obj,hObject)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
%%

labelData = get(gca,'UserData');
labels = labelData.labels;
% frames = labelData.frames;
numColumns = labelData.numColumns;
dims = labelData.dims;
pos = get(event_obj,'Position');

% find the label for each click
% we used row major tiling so flip x<->y
idx = ceil(pos(1)/dims(1));
idy = ceil(pos(2)/dims(1));
idxlabel = sub2ind([numColumns numColumns],idx,idy);

strs = strsplit(labels{idxlabel},' ');
lab = sprintf('idx: %d\nlab: %s\nidxW: %d\nWloc: %.3f,%.3f,%.3f',...
    [idxlabel,...
    strs{1},': ', [strs{3},'-',strs{4},'-',strs{5}]],...
    str2num(strs{7}),...
    str2num(strs{9}),str2num(strs{10}),str2num(strs{11}));

output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
% if length(pos) > 2
output_txt{end+1} = [lab];
disp(output_txt{end})
