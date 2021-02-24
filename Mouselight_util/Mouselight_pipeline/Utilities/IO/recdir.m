function recdir(inputfolder,args,level)
%%
% args.level = opt.level;
% args.ext = opt.ext;
% if exist(opt.seqtemp, 'file') == 2
%     % load file directly
% else
%     args.fid = fopen(opt.seqtemp,'w');
%     recdir(opt.inputfolder,args)
% end
% fid=fopen(opt.seqtemp,'r');
% myfiles = textscan(fid,'%s');
% myfiles = myfiles{1};
% fclose(fid)

% get sequence
if nargin ==0 
    inputfolder = '/nobackup2/mouselight/cluster/SitchedProbability_GN1/2015-06-19-erhan-GN1-maskProb'
    args.level = 5;
    args.fid = fopen('paths.txt','w');
    args.ext = 'tif'
%     try
%         recdir(inputfolder,args,level)
%     catch
%         fclose(args.fid)
%     end
end
%%
if nargin <3
    level = 0;
end
if args.ext(1) == '.'
    args.ext(1) = [];
end
dirinfo = dir(inputfolder);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
tf = cellfun(@(x) x(1)=='.',{dirinfo.name}); %remove hidded and parent folders
% ismember( {dirinfo.name}, {'.', '..','.log','.backup'});
dirinfo(tf) = [];  %remove current and parent directory.

if isfield(args,'pattern')
    matches = regexp({dirinfo(:).name},args.pattern);
    keep_these = zeros(1,length(dirinfo));
    for ifold = 1:length(matches)
        if matches{ifold} & matches{ifold}(1)==1
            keep_these(ifold) = 1;
        end
    end
    dirinfo(~keep_these) = []; %remove any folder that doesn't follow pattern
end
%%
if level == args.level
    % search file
    % get files with argument
    myfiles = dir([inputfolder,'/*.',args.ext]);
    % append to xls file
    for ii=1:length(myfiles)
        if 0&exist(fullfile(inputfolder,[myfiles(ii).name(1:end-3),'.tif']),'file')
            
        else
            fprintf(args.fid,'%s\n',fullfile(inputfolder,myfiles(ii).name));
        end
    end
%     % append to xls file
%     for ii=1:length(myfiles)
%         if myfiles(ii).bytes~=789753938
%             fprintf(args.fid,'%s\n',fullfile(inputfolder,myfiles(ii).name));
%         end
%     end
    %return
else
    % recursion
    for idx=1:length(dirinfo)
        recdir(fullfile(inputfolder,dirinfo(idx).name),args,level+1)
    end
    if level == 1
        disp(sprintf('Finished %d : %s',level,inputfolder))
    end
end
if level==0
    fclose(args.fid);
end

end