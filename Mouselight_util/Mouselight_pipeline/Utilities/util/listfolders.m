function [varargout] = listfolders(varargin)

numfold = length(varargin)/2;
varargout = cell(1,numfold);
for ifold = 1:numfold
    fold = varargin{(ifold-1)*2+1};
    ext = varargin{(ifold)*2}; 
    % list tif files
    clear args
    pathfile = fullfile(fold,sprintf('list%sfiles',ext));

    args.ext = ext;
    args.level = 3;
    if exist(pathfile, 'file') == 2
        % load file directly
    else
        args.fid = fopen(pathfile,'w');
        recdir(fold,args)
        % make it write protected
        unix(sprintf('chmod -w %s',pathfile));
    end
    
    clear args
    fid = fopen(pathfile);
    list = textscan(fid,'%s','Delimiter','\n');
    varargout{ifold}=list{1};
    fclose(fid);
end
end
