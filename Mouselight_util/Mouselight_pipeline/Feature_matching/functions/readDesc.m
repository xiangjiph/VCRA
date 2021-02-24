function des = readDesc(inputfolder,desc_ch,ext_desc)
%%
if nargin<2
    desc_ch = {'1'};
    ext_desc = 'txt';
end
if nargin<3
    ext_desc = 'txt';
end
myfiles = dir(fullfile(inputfolder,['*.',ext_desc]));
fid = fopen(fullfile(myfiles(1).folder,myfiles(1).name));
tLines = fgets(fid);
numch = length(desc_ch);
delimiter = ' '; %or whatever
if length(tLines)>1
    numCols = numel(strfind(tLines,delimiter)) + 1;
    fclose(fid);
else
    des = [];
    return
    % error('EMPTY DESCRIPTOR FILE !!')
end
format = repmat('%f ',1,numCols);format(end)=[];
%% if there are multi-channel descriptor, decide what to do
%% TODO
clear des
thesechannels = desc_ch;

for idx_=1:length(thesechannels)
    idx = thesechannels{idx_};
    ii = str2double(idx)+1;
    myfid1 = fopen(fullfile(inputfolder,myfiles(ii).name));
    data = textscan(myfid1,format); % round locations, there is a bug in estimation that shifts z-locations 0.5 pix. rounding results in better MSE
    fclose(myfid1);
    if size(data,2)==3
        des{ii} = [data{1} data{2} data{3}];
    else
        des{ii} = [[data{1} data{2} data{3}] data{4:numCols}];
    end
end
%%
%% TODO make it more than 2 channel
if numch>1
    pd2 = pdist2(des{1},des{2});
    [aa1,bb1] = min(pd2,[],1);
    [aa2,bb2] = min(pd2',[],1);% instead of min(pd2,[],2), make it row vector with transpose to prevent dimension error for single row entities
    bb1(aa1>1) = 1;
    bb2(aa2>1) = 1;
    keepthese = [1:length(bb1)]==bb2(bb1);
    des = des{2}(keepthese,:);
else
    des = des{end};
end

%%

