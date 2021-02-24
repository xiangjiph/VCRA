function mynames = getRelativePath(mylist)

mynames = cell(1,length(mylist));
parfor ii=1:length(mylist)
    strs = strsplit(mylist{ii},'/');
    mynames{ii} = fullfile(strs{end-3:end});
end
end