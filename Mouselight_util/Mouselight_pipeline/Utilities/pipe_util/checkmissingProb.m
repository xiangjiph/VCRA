function checkmissingProb(inputfolder,outputtfolder)
%DEPLOYONBRAIN Summary of this function goes here
% 
% [OUTPUTARGS] = DEPLOYONBRAIN(INPUTARGS) Explain usage here
% 
% Inputs: 
% 
% Outputs: 
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2016/02/09 09:11:11 $	$Revision: 0.1 $
% Copyright: HHMI 2016
addpath(genpath('./common'))
addpath(genpath('./functions'))
if nargin==2
    logfolder = [];
end
%%
logtag = 'ax-%s-log.%s.txt';
nametag = 'prob'; % output tag, e.g. 0000-<nametag>.0.h5
numcores = 4;
memsize = numcores*7.5*1000;
ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.9-Linux/run_ilastik.sh';
ilpfile = '/groups/mousebrainmicro/mousebrainmicro/erhan_dm11/AxonClassifier/axon_uint16.ilp';
if 0
    outext = 'tif'; %#ok<UNRCH>
    outextformat = '"multipage tiff"';
else
    outext = 'h5';
    outextformat = '"hdf5"';
end

out = outputtfolder;
myshfile = fullfile(outputtfolder,sprintf('cluster_ilastik_%s.sh',date));

mkdir(out);
unix(sprintf('umask g+rxw %s',out));
unix(sprintf('chmod g+rxw %s',out));

if ~isempty(logfolder)
    mkdir(logfolder);
    unix(sprintf('umask g+rxw %s',logfolder));
    unix(sprintf('chmod g+rxw %s',logfolder));
end
%%
% delete h5list as it changes every iteration
h5pathfile = fullfile(out,'listh5files');
unix(sprintf('rm %s',h5pathfile))

clear args
pathfile = fullfile(outputtfolder,'listtiffiles');
args.ext = 'tif';
args.level = 3;
if exist(pathfile, 'file') == 2
    % load file directly
else
    args.fid = fopen(pathfile,'w');
    recdir(inputfolder,args)
    % make it write protected
    unix(sprintf('chmod -w %s',pathfile));
end

% first get file list
filename = pathfile; fid = fopen(filename); targetlist = textscan(fid,'%s','Delimiter','\n');targetlist=targetlist{1}; fclose(fid);
% get the portion after inputfolder
mynames = cell(1,length(targetlist));for ii=1:length(targetlist); mynames{ii} = targetlist{ii}(length(inputfolder)+1:end); end

%% curation
% curationH5(inputfolder,outputtfolder)

%% missing files
if 0
    missingfiles=ones(1,length(targetlist)); %#ok<UNRCH>
elseif 1
    if isempty(logfolder)
        % for every input tile, check the output tile log folder
        missingfiles = checkmissingpertile(out,mynames,nametag,outext,logtag);
    else
        pathfile = fullfile(outputtfolder,'listtiffiles'); %#ok<UNRCH>
        missingfiles = checkmissing(pathfile,logfolder);
    end
end
sum(missingfiles)

%%
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 10;
fid = fopen(myshfile,'w');
for ii=1:1:length(targetlist)
    %%
    if ~missingfiles(ii)
        continue
    end
%     if ~missingfiles(ii) | scopeloc.loc(round(ii/2),2)>44283.7/1e3 | scopeloc.loc(round(ii/2),3)>23474.8/1e3
%         continue
%     end
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    jobname = sprintf('ilp_%05d-%s',ii,randString);
    infiles=targetlist{ii};
    
    [subpath,name,~]=fileparts(targetlist{ii}(length(inputfolder)+1:end));
    % TODO: rename using prob 
    subname = strsplit(name,{'-','.'});
    if ~exist(fullfile(out,subpath), 'dir')
        mkdir(fullfile(out,subpath));
        unix(sprintf('chmod g+rwx %s',fullfile(out,subpath)));
    end
    if ~exist(fullfile(out,subpath,'.log'), 'dir')
        mkdir(fullfile(out,subpath,'.log'));
        unix(sprintf('chmod g+rwx %s',fullfile(out,subpath,'.log')));
    end
    outputformat=fullfile(out,subpath,sprintf('%s-%s.%s.%s',subname{1},nametag,subname{3},outext));
    logfile = fullfile(out,subpath,'.log',sprintf(logtag,subname{1},subname{3}));
    if ~isempty(logfolder)
        logfile=fullfile(logfolder,sprintf('ilp_%05d-%s.txt',ii,randString));
    end
    
    argsout = sprintf('''%s --headless  --cutout_subregion="[(None,None,None,0),(None,None,None,1)]" --logfile=%s --project=%s --output_format=%s --output_filename_format=%s %s''',...
        ilastikloc,logfile,ilpfile,outextformat,outputformat,infiles);
    mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d bsub -n%d -We 25 -J t-%d-%s -o /dev/null %s\n',numcores,memsize,numcores,ii,jobname,argsout);
    fwrite(fid,mysub);
end
fclose(fid);
unix(sprintf('chmod +x %s',myshfile));


end

function missingfiles = checkmissingpertile(out,mynames,nametag,outext,logtag)
%%
missingfiles=ones(1,length(mynames));
parfor idx = 1:length(mynames)
    outputfile = fullfile(out,strrep(strrep(mynames{idx},'ngc',nametag),'tif',outext));
    [mypath,myfile,myext] = fileparts(outputfile);
    mychannel = myfile(end);
    logfile = fullfile(mypath,'.log',sprintf(logtag,myfile(1:5),mychannel));
    if exist(outputfile,'file')
        % check if log file has completed line
        if exist(logfile,'file')
            [q,w] = system(['grep "Completed Batch Processing" ',logfile]);
            missingfiles(idx) = q;
        end
    end
end
end

function missingfiles = checkmissing(pathfile,logfolder)
%%
filename = pathfile;
fid = fopen(filename);
targetlist = textscan(fid,'%s','Delimiter','\n');targetlist=targetlist{1};
fclose(fid);
myfold = logfolder;
myfiles = dir(sprintf('%s/*.txt',myfold))
%%
valind = zeros(1,length(targetlist));
valfiles = zeros(length(myfiles),2);
parfor_progress(length(myfiles));
parfor ii=1:length(myfiles)
    parfor_progress;
    if ~rem(ii,1e3)
        ii
    end
    myfile = fullfile(myfold,myfiles(ii).name);
    [~,mf,~] = fileparts(myfile);
    [aa] = strsplit(mf,'_');
    tar = str2num(aa{2}(1:5));
    %%
    %%
    [q,w] = system(['tail -n ',num2str(1),' ',myfile]);
    if length(w)<27
        valfiles(ii,:) = [tar, 0];
    else
        valfiles(ii,:)=[tar strcmp(w(end-26:end-1),'Completed Batch Processing')];
    end
%     if length(w)<27
%         valind(tar)=0;
%     else
%         valind(tar)=valind(tar)+strcmp(w(end-26:end-1),'Completed Batch Processing');
%     end
end
parfor_progress(0);
%%
missingfiles = zeros(1,length(targetlist));
missingfiles(setdiff(1:length(targetlist),valfiles(valfiles(:,2)>0,1))) = 1;
%%


end
function recdir(inputfolder,args,level)
%%
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
tf = ismember( {dirinfo.name}, {'.', '..'});
dirinfo(tf) = [];  %remove current and parent directory.
%%
if level == args.level
    % search file
    % get files with argument
    myfiles = dir([inputfolder,'/*.',args.ext]);
    % append to xls file
    for ii=1:length(myfiles)
        fprintf(args.fid,'%s\n',fullfile(inputfolder,myfiles(ii).name));
    end
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

function missmatch()
inputfolder = '/nrs/mouselight/cluster/classifierOutputs/2015-06-19_backup/'
filename = fullfile(inputfolder,'listfiles');
fid = fopen(filename);
targetlist = textscan(fid,'%s','Delimiter','\n');targetlist=targetlist{1};
fclose(fid);
myfun = @(x) strsplit(x,'/');
clear mynames
for ii=1:length(targetlist)
    %     K = myfun(targetlist{ii});
    % get the portion after inputfolder
    mynames{ii} = targetlist{ii}(length('/groups/mousebrainmicro/mousebrainmicro/from_tier2/data/2015-06-19/Tiling')+1:end);
end
%%
inputfolder = '/nrs/mouselight/cluster/classifierOutputs/2015-06-19_backup/classifier_output'
filename = fullfile(inputfolder,'filelist.txt');
fid = fopen(filename);
C2 = textscan(fid,'%s','Delimiter','\n');C2=C2{1};
fclose(fid);
myfun = @(x) strsplit(x,'/');
clear mynames2
for ii=1:length(C2)
    %     K = myfun(targetlist{ii});
    % get the portion after inputfolder
    mynames2{ii} = C2{ii}(length(inputfolder)+1:end);
end
%%
found=zeros(1,size(mynames,2));
iis1 = [2:27,32];
iis2 = [2:27,33];
for ii=1:size(mynames,2)
    %%
    ii
    % check if this exists
    mystr = mynames{ii}(iis1);
    for jj=1:size(mynames2,2)
        mystr2=mynames2{jj}(iis2);
        if all(mystr==mystr2)
            found(ii)=1;
            mynames2(jj)=[];
            break
        end
    end
end
end




