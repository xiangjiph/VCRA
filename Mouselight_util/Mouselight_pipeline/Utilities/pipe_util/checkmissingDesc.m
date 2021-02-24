function checkmissingDesc(descinput,descoutput)
if nargin ==1
    descoutput = descinput;
end
args.level = 3;
args.ext = 'h5';

opt.inputfolder = descinput;
opt.seqtemp = fullfile(opt.inputfolder,'listh5files')
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
inputfiles = textscan(fid,'%s');
inputfiles = inputfiles{1};
fclose(fid);
if 1
    missingfiles=ones(1,length(inputfiles));
    % for every tif check if h5 exists
    parfor ii=1:length(inputfiles)
        if ~missingfiles(ii)
            continue
        end
        if exist(strrep(strrep(inputfiles{ii},args.ext,'txt'),'prob','desc'),'file')
            missingfiles(ii) = 0;
        end
    end
else
    missingfiles = ones(1,size(inputfiles,1));
end
if ~sum(missingfiles)
    disp('Found match for all tiles')
    return
else
    sum(missingfiles)
end
    
%% mcc -m -R -nojvm -v <function.m> -d <outfolder/>  -a <addfolder>
numcores = 4;
imagesiz = [1024 1536 251];
pre = 'prob' % or ngc
post = 'desc'
rt = 4;
myfile = fullfile(pwd,sprintf('./shfiles/dogdescriptorrun_%s_miss.sh',date));
% myshfile = fullfile(experimentfolder,sprintf('cluster_ilastik_%s.sh',brain));

compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/dogDescriptor/dogDescriptor'
if 0
    mkdir(fileparts(compiledfunc))
    unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
    sprintf('mcc -m -v -R -singleCompThread %s/dogDescriptor.m -d %s',pwd,fileparts(compiledfunc))
    unix(sprintf('chmod g+rwx %s',compiledfunc))
end

%find number of random characters to choose from
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
numRands = length(s);
%specify length of random string to generate
sLength = 10;
ROI = [5 imagesiz(1)-5 5 imagesiz(2)-5 5 imagesiz(3)-1];
%-o /dev/null
esttime = 6*60;
%%
fid = fopen(myfile,'w');
% outputlogfold = fullfile('/groups/mousebrainmicro/mousebrainmicro/LOG',brain,'descriptor')
% mkdir(outputlogfold)
% unix(sprintf('chmod g+rwx %s',outputlogfold))
logout=0
vv=find(missingfiles);
% vv(vv<25497)
for ii=vv;%(vv<22131)%find(missingfiles)%
    %%
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('dog_%05d-%s',ii,randString);
    if 1
        [file_path,file_name,file_ext] = fileparts(inputfiles{ii});
        outfile = fullfile(file_path,[strrep(file_name,pre,post),'.txt']);
    else
        outfile = fullfile(outputfold,sprintf('%05d-%s.%d.txt',floor((ii-1)/2)+1,pre,rem(ii+1,2)));
        %     outfile = fullfile(outputfold,sprintf('%05d-%s.%d.txt',((ii)),pre,0));
    end
    if 0
        if logout
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d> %s/output-%05d.log''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt,outputlogfold,ii);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,outputlogfold,argsout);
        else
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,'/dev/null',argsout);
        end
    else
        if logout
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d> %s/output-%05d.log''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt,outputlogfold,ii);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o %s -b y -cwd -V %s\n',numcores,esttime,name,outputlogfold,argsout);
        else
            argsout = sprintf('''%s %s %s "[%d %d %d]" "[%f %f %f]" "[%f %f %f]" "[%d %d %d %d %d %d]" %d''',compiledfunc,inputfiles{ii},outfile,...
                11*ones(1,3),3.4055002*ones(1,3),4.0498447*ones(1,3),ROI,rt);
            mysub = sprintf('bsub -n%d -We %d -J %s -o %s %s\n',numcores,esttime/60,name,'/dev/null',argsout);
        end
    end
%     dogDescriptor(inputfiles{ii},outfile,sprintf('[%d %d %d]',11*ones(1,3)),sprintf('[%f %f %f]',3.4055002*ones(1,3)),sprintf('[%f %f %f]',4.0498447*ones(1,3)),...
%         sprintf('[%d %d %d %d %d %d]',ROI),sprintf('%d',rt))
    fwrite(fid,mysub);
end
unix(sprintf('chmod +x %s',myfile));
fclose(fid);
sprintf('%s',myfile)