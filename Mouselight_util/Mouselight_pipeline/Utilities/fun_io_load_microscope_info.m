function params = fun_io_load_microscope_info(acqusitionfolder)

myfiles = dir(fullfile(acqusitionfolder,'*.acquisition'));
fid = fopen(fullfile(myfiles(1).folder,myfiles(1).name));
delimiter = ':'; %or whatever
while true
    tLines = fgetl(fid);
    if tLines<0
        break
    end
    strs = strsplit(tLines,delimiter);
    if length(strs)>1
        params.(strip(strs{1})) = str2double(strs{2});
    end
end
fclose(fid);
end