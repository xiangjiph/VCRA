function params = readScopeFile(acqusitionfolder)
myfiles = dir(fullfile(acqusitionfolder,['*.acquisition']));
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

myfiles = dir(fullfile(acqusitionfolder,['*.microscope']));
fid = fopen(fullfile(myfiles(1).folder,myfiles(1).name));
delimiter = ':'; %or whatever
while true
    tLines = fgetl(fid);
    if tLines<0
        break
    elseif length(tLines)>2 & strcmp(tLines(1:3),'fov')
        while ~strcmp(tLines,'}')
            tLines = fgetl(fid);
            strs = strsplit(tLines,delimiter);
            if length(strs)>1
                params.(strip(strs{1})) = str2double(strs{2});
            end
        end
        break
    end
end
fclose(fid);