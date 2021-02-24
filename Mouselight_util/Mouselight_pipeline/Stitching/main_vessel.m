function [outputArgs] = main_vessel(brain,inputfolder,pipelineoutputfolder)
%STICHING pipeline. Reads scope generated json file and returns a yml
%configuration file that goes into renderer. Requires calling cluster jobs
%to create subresults, i.e. descriptors. These functions can also run in
%local machines with proper settings.
%
% [OUTPUTARGS] = STICHING(jsonfile)
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

% $Author: base $	$Date: 2016/09/21 11:52:40 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%% Pipeline test
brain = '2020-02-01';
% inputfolder = '/nrs/mouselight/pipeline_output/2020-02-01/stage_1_line_fix_output';
inputfolder = fullfile('/groups/mousebrainmicro/mousebrainmicro/data', brain, 'Tiling');
pipelineoutputfolder = '/nrs/mouselight/pipeline_output/2020-02-01';
runfull = true;
%% MAKE SURE PATHS etc are correct
% inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain);
runfull = false;
if nargin==1
    brain = '2018-08-01';
    inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',brain);
    pipelineoutputfolder = sprintf('/nrs/mouselight/pipeline_output/%s',brain)
elseif nargin<1
    error('At least pass brain id')
end
%% Set directories
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
directions = 'Z';

addpath(genpath('./common'))
addpath(genpath('./functions'))
tag='';
% classifierinput = inputfolder;
% raw input to descriptor generotion

piperun = 1;

if piperun
    %         classifieroutput = fullfile(pipelineoutputfolder,'stage_2_classifier_output')
    %         descinput = classifieroutput;
    descoutput = fullfile(pipelineoutputfolder,'stage_2_descriptor_output');
    matchinput = descoutput;
    matchoutput = fullfile(pipelineoutputfolder,'stage_3_point_match_output');
end

experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);
% experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s_xj_downsample',brain);
% experimentfolder = '/nrs/mouselight/cluster/classifierOutputs/2018-08-15_xj_wholebrain';
matfolder = fullfile(experimentfolder,'matfiles/');
if ~isfolder(matfolder)
    mkdir(matfolder)
    unix(sprintf('umask g+rxw %s',matfolder));
    unix(sprintf('chmod g+rxw %s',matfolder));
end
scopefile = fullfile(matfolder,'scopeloc.mat');
if piperun
    descriptorfolder = descoutput;
    matchfolder = matchoutput;
else
    descriptorfolder = fullfile(experimentfolder,'classifier_output');
    matchfolder = descriptorfolder;
end

desc_ch = {'0'};
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file
matchedfeatfile = fullfile(matfolder,sprintf('feats_ch%s.mat',desc_ch{:})); % accumulated descriptor file
%% 0: INTIALIZE
% read scope files, determine the valid neighbor indices in the file list
% and save the information to scopefile
if runfull
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    % scopeloc.gridix = [x y z cut_count]; grid subscripts
    % scopeloc.loc: [x, y, z]: stage position in mm
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); 
    % neighbors: N-by-7 array, [idx -x -y +x +y -z +z] indices of neighbors
    % of tile(idx)
    save(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
end
% %% BULLSHIT CURATION STUFF
% % Are we going to use this part again? 
% % obsolute after pipeline, TODO: fix missing condition for tile runs
% % rather then channel logs
% if 0
%     curationH5(classifierinput,classifieroutput)
% end
% 
% if 0
%     % checkmissingProb(classifierinput,classifieroutput)
%     checkmissingDesc(descinput,descoutput)
%     checkmissingMatch(matchinput,matchoutput)
% end
%% LOAD MATCHED FEATS
if runfull
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder');
    directions = 'Z';
    checkversion = 1; % 1: loads the version with "checkversion" extension and overwrites existing match if there are more matched points
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,matchfolder,directions,checkversion);
    % regpts: structure array with fields from descriptor match and
    % neighboring block index list
    % featmap: structure array with fields: scopefile 1, scopefile 2,
    % paireddescriptor 
    save(fullfile(matfolder,'regpts.mat'),'regpts','featmap', '-v7.3')
    if ~exist(fullfile(matfolder,'regpts_1stiter.mat'),'file') % faster to make a copy
        unix(sprintf('cp %s %s',fullfile(matfolder,'regpts.mat'),fullfile(matfolder,'regpts_1stiter.mat')))
    end
end
% %% Re-run the feature matching algorithm
% % For improving the feature matching 
% if 0 % iterate on missing tiles (ANOTHER BULLSHIT)
%     
%     addpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'),'-end')
%     %pointmatch_task(brain,runlocal)
%     directions = 'Z';
%     ch=desc_ch{1};
%     [~,sample] = fileparts(experimentfolder);
%     runlocal=1;
%     pointmatch_task_local_vessel(sample,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal)
%     rmpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'))
% end
%% Process the matched feature - Main part of stitching algorithm
% 2 scope params estimation.
% i) finds matches on x&y
% ii) finds field curvature based on matched points
% iii) creates a 3D affine model by jointly solving a linear system of
% equations
processing_log_fp = sprintf('~/stitching_log_2019_05_07.txt');
if runfull | 1
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    % paramater setting for descrtiptor match
    scopeacqparams = readScopeFile(fileparts(scopeloc.filepath{1}));
    % params.sample = brain;
    params.scopeacqparams = scopeacqparams;
    params.viz = 0;
    params.debug = 0;
    params.Ndivs = 4;
    params.Nlayer = 4;
    params.htop = 5;
    params.expensionratio = 1;
    params.imagesize = [1024 1536 251];
    params.imsize_um = [scopeacqparams.x_size_um scopeacqparams.y_size_um scopeacqparams.z_size_um];
    params.overlap_um = [scopeacqparams.x_overlap_um scopeacqparams.y_overlap_um scopeacqparams.z_overlap_um];
    params.order = 1;
    params.applyFC = 1;
    params.beadparams = [];%PLACEHOLDER FOR BEADS, very unlikely to have it...
    params.singleTile = 1;
    if 1
        diary(processing_log_fp);
        diary('on');
        [paireddescriptor,curvemodel,scopeparams] = ...
            tileProcessor_vessel(scopeloc,descriptorfolder,desc_ch,params);
        %             save(descriptorfile,'descriptors','-v7.3')
        save(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
            'scopeparams', 'curvemodel','params','-v7.3')
        diary('off');
    else
        descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);
        [paireddescriptor, scopeparams, R, curvemodel,scopeparams_, paireddescriptor_,curvemodel_] = ...
            estimateFCpertile(descriptors,neighbors,scopeloc,params);
        save(descriptorfile,'descriptors','-v7.3')
        save(fullfile(matfolder,'scopeparams_pertile_singletile'),'paireddescriptor', ...
            'scopeparams', 'R', 'curvemodel','scopeparams_', 'pareddescriptor_', ...
            'curvemodel_','params')
    end
end

%% Generate descriptor match quality heat map
if runfull
    %%
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'scopeparams_pertile'),'scopeparams')
    load(fullfile(matfolder,'regpts'),'regpts')
    output_fig_root = fullfile('~/Documents/GitHub/stitching_quality_vis/', ...
        sprintf('%s_%s', brain, date));
%     videofile = sprintf('./videos/%s-1stiter-ch1-%s',brain,date)
    % descriptorMatchQuality(regpts,scopeparams,scopeloc,videofile)
    %     createThumb(regpts,scopeparams,scopeloc,videofile)
    visualization_map_on_the_fly = false;
    descriptorMatchQualityHeatMap_vessel(regpts,scopeparams{end},scopeloc,output_fig_root)
%     descriptorMatchQualityHeatMap_forPaper(regpts,scopeparams{end},scopeloc,videofile)
end
%% Compute the vector field
tic
if 1
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'regpts'),'regpts')
    load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
        'scopeparams', 'curvemodel','params')
    
    vecfield3D = vectorField3D(params,scopeloc,regpts,scopeparams{end},curvemodel{end},[]);
    if 1
        save(fullfile(matfolder,sprintf('%s_%s',datestr(now,'mmddyyHHMMSS'),'vecfield3D')),'vecfield3D','params')
        save(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
    end
end
toc
load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
load(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
vecfield = vecfield3D;

% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
params.big = 1;
params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
sub = 0;
params.root = vecfield.root;

if sub
    targetidx = getTargetIDx(scopeloc,neighbors);
    params.outfile = fullfile(experimentfolder,sprintf('%s_sub.control.yml',date));
else
    targetidx = 1:size(scopeloc.gridix,1);
    params.outfile = fullfile(experimentfolder,sprintf('%s.control.yml',date));
end

writeYML(params, targetidx(:)', vecfield)
unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache.yml')))
%
if ~sub
    params.big=0
    params.outfile = sprintf('%s/%s.old.control.yml',experimentfolder,date);
    writeYML(params, targetidx(:)', vecfield)
    unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache_old.yml')))
end
return
%% 0.1: FLAT RUN
% generate yml for without any optimization
if 0
    %%
    load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
    if scope==1
        scope1_beadparams = load('./beadparams/scope1_beadparams');
        scopeparams = scope1_beadparams.scope1_beadparams;
    else
        scope2_beadparams = load('./beadparams/scope2_beadparams');
        scopeparams = scope2_beadparams.scope2_beadparams;
    end
    vecfield = vectorField_flatrun(params,scopeloc,scopeparams,2);
    
    load ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
    [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
        beadparams,neighbors,scopeloc,paireddescriptor,R,curvemodel,imsize_um);
    vecfield3D_flat_4neig = vectorField_flatrun_pertile(params,scopeloc,scopeparams_,curvemodel_,[]);
    save pertile_4neig scopeparams scopeparams_ paireddescriptor_ curvemodel_ vecfield3D_flat_4neig
end
%% stitching quality test
if 0
    load(fullfile(matfolder,'scopeloc'),'scopeloc','imsize_um','experimentfolder','inputfolder')
    load(fullfile(matfolder,'vecfield'),'vecfield','params')
    %%
    clc
    params.big = 1;
    params.dims = [params.imagesize 2]%[1024 1536 251 2]
    sub = 0;
    inds_ = inds(1)';
    neigs = neighbors(inds_,checkthese);
    targetidx = neigs([1 3])
    params.root = vecfield.root;
    %%
    %
    params.outfile = sprintf('%s%s-%d_%d_sub_1.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(1)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y1-ccx2-fixxed')
    %%
    
    params.outfile = sprintf('%s%s-%d_%d_sub_2.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(2)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y2-ccx2-fixxed')
end
end
