% compares h5 and raw folders and mv extra folders in h5 into backup
% folder. run this first, then check for missing tiles
function curationH5(tiffolder_,h5folder_)
% sample = '2017-09-25'
% logfolder = '/groups/mousebrainmicro/mousebrainmicro/LOG/pipeline/'
% axlogfiles = dir(fullfile(logfolder,['ax-2017-08*.txt']));
% h5folder='/nrs/mouselight/cluster/classifierOutputs/';
% tiffolder='/groups/mousebrainmicro/mousebrainmicro/data';
% h5folder_ = fullfile(h5folder,sprintf('%s/classifier_output',sample));
% tiffolder_ = fullfile(tiffolder,sprintf('%s/Tiling',sample));

% move missing folders into backup
backupfolder = fullfile(h5folder_,'.backup');

[inputtiflist,inputh5list] = listfolders(tiffolder_,'tif',h5folder_,'h5');

%%
tifpaths = getRelativePath(inputtiflist);
h5paths = getRelativePath(inputh5list);
%% 
% compare paths
numh5 = length(h5paths);
numtif = length(tifpaths);

idxh5 = zeros(1,numh5);
parfor ii=1:numh5
    % search in tif list
    found = find(contains(tifpaths,strrep(strrep(h5paths{ii},'prob','ngc'),'h5','tif')));
    if found
        idxh5(ii) = found;
    end
end
% folder list
folderpaths2backup = unique(cellfun(@fileparts,{inputh5list{idxh5==0}},'UniformOutput',false));
disp(folderpaths2backup)
%%
reply = input(sprintf('Found %d mismatching folders, do you want to curate ? Y/N [Y]:',length(folderpaths2backup)),'s');
if isempty(reply)
    reply = 'N';
end
if reply== 'Y' | reply=='y'
    for ii=1:length(folderpaths2backup)
        strs = strsplit(folderpaths2backup{ii},'/');
        mkdir(fullfile(backupfolder,fullfile(strs{end-2:end})))
        unix(sprintf('mv %s %s',folderpaths2backup{ii},fullfile(backupfolder,fullfile(strs{end-2:end-1}))));
    end
    % update h5 filelist
    unix(sprintf('rm %s',fullfile(h5folder_,'listh5files')));
    [inputh5list] = listfolders(h5folder_,'h5');
else
    disp('Skipping curation')
end

