function [outputArgs] = descriptorMatchQualityHeatMap_vessel(regpts,scopeparams,scopeloc,output_folder)
%DESCRIPTORMATCHQUALITY Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCHQUALITY(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/11/15 11:13:34 $	$Revision: 0.1 $
% Copyright: HHMI 2016

%%
visualize_Q = false;
% populate afftile
numTiles = size(regpts,2);
afftile = zeros(3, 4, numTiles);
medaffine = zeros(3,3,numTiles);
for ii=1:numTiles
    if isempty(scopeparams(ii).affineglFC)
        continue;
    end
    medaffine(:,:,ii) = scopeparams(ii).affineglFC;
end
afftransform = median(medaffine,3);
mflip = -eye(3);%mflip(3,3)=1;
afftransform = afftransform * mflip;
for idxt = 1:numTiles
    %%
    % form an affine matrix
    tform = [afftransform scopeloc.loc(idxt,:)'*1e6];
    afftile(:,:,idxt) = tform;
end
%%
% D is the distance between two z-adjacent tiles recorded by the stage
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
checkthese = [1 4 5 7]; % 0 - below
D = zeros(1,size(neighbors,1));
for ii=1:size(neighbors,1)
    if any(isnan(neighbors(ii,[1 7])))
        D(ii) = nan;
        continue
    end
    D(ii)=diff(scopeloc.loc(neighbors(ii,[1 7]),3));
end

%%
latticeZRange = unique(scopeloc.gridix(:,3));
clear xyz_t xyz_tp1 XYZ_t XYZ_tp1 ranXY ranXYZ
for t = latticeZRange(1:end-1)'
    %%
    disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
    ix = (scopeloc.gridix(:,3)'==t);
    %%
    cnt = 0;
    Npts = 0;
    clear layer layerp1 poslayer_t poslayer_tp1
    [layer layerp1 poslayer_t poslayer_tp1]=deal([]);
    
    for id_ix = find(ix)
        %%
        if isempty(regpts{id_ix}.X) | isnan(regpts{id_ix}.neigs(end))
            vec{t} = [];
            continue
        end
        %%
        cnt = cnt+1;
        bb(cnt)=id_ix;
        neigs = regpts{id_ix}.neigs;
        layer{cnt} = regpts{id_ix}.X; % Pixel position in image stack
        layerp1{cnt} = regpts{id_ix}.Y;
        % Layer_t and Lyaer_tp1 are the position of matched descriptors
        % after affine transformation for correcting the homography. 
        Layer_t = [regpts{id_ix}.X, ones(size(regpts{id_ix}.X,1),1)] * afftile(:,:,neigs(1))';
        Layer_tp1 = [regpts{id_ix}.Y, ones(size(regpts{id_ix}.Y,1),1)] * afftile(:,:,neigs(4))';
        Npts = Npts + size(Layer_t,1);
        poslayer_t{cnt} = Layer_t; % Physical position
        poslayer_tp1{cnt} = Layer_tp1;
        dL = Layer_t - Layer_tp1; % XY diplacement between adjacent tiles
        ranXY(t) = max(sqrt(sum(dL(:,1:2).^2, 2))) / 1e3; % Euclidean distance between matched descriptors in z-adjecent tiles on xy plane
        ranXYZ(t) = max(sqrt(sum(dL(:,1:3).^2, 2))) / 1e3; % Euclidean distance between matched descriptors in z-adjecent tiles in 3D splace
        vec{t}{id_ix} = median(dL(:,1:3));
    end
    %%
    if isempty(layer)
        continue
    end
    xyz_t{t} = cat(1,layer{:});
    xyz_tp1{t} = cat(1,layerp1{:});
    XYZ_t{t} = cat(1,poslayer_t{:});
    XYZ_tp1{t} = cat(1,poslayer_tp1{:});
end
%%
centerPoint = 25;
scalingIntensity = 4;
dataMax = max(ranXY(ranXY>0));
dataMin = min(ranXY(ranXY>0));
newMap = newColMap(centerPoint,scalingIntensity,dataMin,dataMax);

trmax=cellfun(@max, XYZ_t(latticeZRange(1):end), 'UniformOutput', false);
trmin=cellfun(@min, XYZ_t(latticeZRange(1):end), 'UniformOutput', false);
Rmax = max(cat(1,trmax{:})) + scopeparams(1).imsize_um * 1e3 .* [1 1 0];
Rmin = min(cat(1,trmin{:})) - scopeparams(1).imsize_um * 1e3 .* [1 1 0];
% plot desctiptors
% myfig = 100;
% figure(myfig), cla, clf
% hold on
loops = latticeZRange(end)-latticeZRange(1);
F(loops) = struct('cdata',[],'colormap',[]);
iter=1;
%%
for t = latticeZRange(1:end-1)'%[779,780]%
    %%
    ix = (scopeloc.gridix(:,3)'==t);
    x = scopeloc.loc(ix,1)*1e6-scopeparams(1).imsize_um(1)*1e3;
    y = scopeloc.loc(ix,2)*1e6-scopeparams(1).imsize_um(1)*1e3;
    w = scopeparams(1).imsize_um(1)*ones(sum(ix),1)*1e3;
    h = scopeparams(1).imsize_um(2)*ones(sum(ix),1)*1e3;
    clear theseinds
    theseinds = 1:sum(ix);
    if t>size(XYZ_t,2)
        F(iter) = getframe;
        iter=iter+1;
        continue
    end
    emptyslice = isempty(XYZ_t{t});
    %%
    if emptyslice
        [X,Y,Z,xyz,valsXY,valsXYZ,U,V,W] = deal([]);
    else
        X = XYZ_t{t}(:,1);
        Y = XYZ_t{t}(:,2);
        Z = XYZ_t{t}(:,3);
        U = XYZ_tp1{t}(:,1)-XYZ_t{t}(:,1);
        V = XYZ_tp1{t}(:,2)-XYZ_t{t}(:,2);
        W = XYZ_tp1{t}(:,3)-XYZ_t{t}(:,3);
        
        MXYZ = sqrt(U.^2 + V.^2 + W.^2)/1e3;
        MXY = sqrt(U.^2+V.^2)/1e3;
        % Interpolated xy displacement of matched descriptor after affine
        % transformation on the slice
        % Basically this is the nonlinear local displacement due to
        % sectioning
        FxXY = scatteredInterpolant(XYZ_tp1{t},MXY,'linear','nearest');
        % Interpolated xyz displacement on the slice
        FxXYZ = scatteredInterpolant(XYZ_tp1{t},MXYZ,'linear','nearest');
        xyz = scopeloc.loc(ix,:)*1e6;
        valsXY = FxXY(xyz);
        valsXYZ = FxXYZ(xyz);
    end
    %% initalize canvas
    if visualize_Q
        fig_handle = figure('Position', [1024, 1024, 512, 1024], 'Visible', 'on');
    else
        fig_handle = figure('Position', [1024, 1024, 512, 1024], 'Visible', 'off');
    end
    h1=subaxis(2, 1, 1, 'sh', 0.03, 'sv', -0.03, 'padding', .04, 'margin', 0);
    cla, hold on
    disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
    caxis([min(ranXYZ(ranXYZ>0)) max(ranXYZ)*.9])
    set(h1, 'XTick', []);
    set(h1, 'YTick', []);
    set(h1,'Color',[1 1 1]*.5)
    set(h1,'Ydir','reverse')
    set(h1,'Box','on')
    set(h1, ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 10        );
    %     text(Rmin(1)+5e5,Rmin(2)+7e5,'Z',...
    %         'FontSize',40,'Color','k','HorizontalAlignment','left')
    view(0,90)
    colorbar
    colormap(newMap)
    title('Z displacement');
    for ii=theseinds
        rectangle('Position', [x(ii) y(ii) w(ii) h(ii)],'EdgeColor','r')
        if ~emptyslice
            xp = [x(ii) x(ii) x(ii)+w(ii) x(ii)+w(ii)];
            yp = [y(ii) y(ii)+h(ii) y(ii)+h(ii) y(ii)];
            patch(xp,yp,sqrt(valsXYZ(ii)^2-valsXY(ii)^2))
        end
    end
    if ~emptyslice
        scatter(XYZ_t{t}(:,1),XYZ_t{t}(:,2),2,'filled', ...
            'MarkerFaceAlpha',.2,'MarkerFaceColor',[1 1 1]*.5)
    end
    axis equal
    xlim([Rmin(1) Rmax(1)])
    ylim([Rmin(2) Rmax(2)])
    %%
    h2=subaxis(2, 1, 2, 'sh', 0.03, 'sv', -0.03, 'padding', 0.04, 'margin', 0);
    cla, hold on
    set(h2, 'XTick', []);
    set(h2, 'YTick', []);
    set(h2,'Color',[1 1 1]*.5)
    set(h2,'Ydir','reverse')
    set(h2,'Box','on')
    set(h2, ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 10        );
    view(0,90)
    axis equal
    colorbar
    caxis([min(ranXY(ranXY>0)) max(ranXY)*.9])
    title('XY dispacement');
    for ii=1:sum(ix)
        rectangle('Position', [x(ii) y(ii) w(ii) h(ii)])
        if emptyslice
            continue
        end
        
        xp = [x(ii) x(ii) x(ii)+w(ii) x(ii)+w(ii)];
        yp = [y(ii) y(ii)+h(ii) y(ii)+h(ii) y(ii)];
        patch(xp,yp,valsXY(ii))
        %         if length(vec{t})>=fix(ii)
        %             vv=vec{t}{fix(ii)};
        %             if ~isempty(vv)
        %                 vv=normr(vv(1:2));
        %                 %                 quiver(x(ii)+w(ii)/2,y(ii)+h(ii)/2,vv(1),vv(2),w(ii)/4,...
        %                 %                     'LineWidth',2,'Color','m','MaxHeadSize',4)
        %
        %                 headWidth = 6;
        %                 headLength = 3;
        %                 ah = annotation('arrow',...
        %                     'headStyle','plain',...
        %                     'HeadLength',headLength,'HeadWidth',headWidth,...
        %                     'Color','r');
        %                 set(ah,'parent',gca);
        %
        %                 x0 = x(ii)+w(ii)/2;
        %                 y0 = y(ii)+h(ii)/2;
        %                 vx0 = vv(1)*w(ii)/3;
        %                 vy0 = vv(2)*h(ii)/3;
        %
        %                 xn0 = (x0-Rmin(1))/(Rmax(1)-Rmin(1));
        %                 yn0 = (y0-Rmin(2))/(Rmax(2)-Rmin(2));
        %                 vnx0 = vx0/(Rmax(1)-Rmin(1));
        %                 vny0 = vy0/(Rmax(2)-Rmin(2));
        %                 set(ah,'position',[x(ii)+w(ii)/2,y(ii)+h(ii)/2,vv(1)*w(ii)/3,vv(2)*h(ii)/3]);
        %                 set(ah,'position',[xn0,yn0,vnx0,vny0]);
        %             end
        %         end
    end
    if ~emptyslice
        scatter(XYZ_t{t}(:,1),XYZ_t{t}(:,2),2,'filled', ...
            'MarkerFaceAlpha',.2,'MarkerFaceColor',[1 1 1]*.5)
    end
    
    %     text(Rmin(1)+5e5,Rmin(2)+7e5,'Lateral',...
    %         'FontSize',40,'Color','k','HorizontalAlignment','left')
    text(Rmin(1)+1e5,Rmax(2)-7e5,num2str(t),...
        'FontSize',40,'Color','k','HorizontalAlignment','left')
    
    xlim([Rmin(1) Rmax(1)])
    ylim([Rmin(2) Rmax(2)])
    drawnow
    F(iter) = getframe(fig_handle);
    iter=iter+1;
    fig_name = sprintf('section_%d_gridz_%d.png', iter, t);
    fig_fp = fullfile(output_folder, fig_name);
    fig_folder = fileparts(fig_fp);
    if ~isfolder(fig_folder)
        mkdir(fig_folder);
    end
    print(fig_handle, '-dpng', '-r600', fig_fp);
    
end
%%
fprintf('Write vedio/n');
v = VideoWriter(sprintf('%s/video', output_folder),'Motion JPEG AVI');
% v.CompressionRatio = 3;
v.Quality = 100;
v.FrameRate=2;
open(v)
writeVideo(v,F)
close(v)
%%
end

function h=subaxis(varargin)
%SUBAXIS Create axes in tiled positions. (just like subplot)
%   Usage:
%      h=subaxis(rows,cols,cellno[,settings])
%      h=subaxis(rows,cols,cellx,celly[,settings])
%      h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])
%
% SETTINGS: Spacing,SpacingHoriz,SpacingVert
%           Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
%           Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
%           Holdaxis
%
%           all units are relative (i.e. from 0 to 1)
%
%           Abbreviations of parameters can be used.. (Eg MR instead of MarginRight)
%           (holdaxis means that it wont delete any axes below.)
%
%
% Example:
%
%   >> subaxis(2,1,1,'SpacingVert',0,'MR',0);
%   >> imagesc(magic(3))
%   >> subaxis(2,'p',.02);
%   >> imagesc(magic(4))
%
% 2001-2014 / Aslak Grinsted  (Feel free to modify this code.)

f=gcf;



UserDataArgsOK=0;
Args=get(f,'UserData');
if isstruct(Args)
    UserDataArgsOK=isfield(Args,'SpacingHorizontal')&isfield(Args,'Holdaxis')&isfield(Args,'rows')&isfield(Args,'cols');
end
OKToStoreArgs=isempty(Args)|UserDataArgsOK;

if isempty(Args)&&(~UserDataArgsOK)
    Args=struct('Holdaxis',0, ...
        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
        'rows',[],'cols',[]);
end
Args=parseArgs(varargin,Args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});

if (length(Args.NumericArguments)>2)
    Args.rows=Args.NumericArguments{1};
    Args.cols=Args.NumericArguments{2};
    %remove these 2 numerical arguments
    Args.NumericArguments={Args.NumericArguments{3:end}};
end

if OKToStoreArgs
    set(f,'UserData',Args);
end


switch length(Args.NumericArguments)
    case 0
        return % no arguments but rows/cols....
    case 1
        if numel(Args.NumericArguments{1}) > 1 % restore subplot(m,n,[x y]) behaviour
            [x1 y1] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(1)); % subplot and ind2sub count differently (column instead of row first) --> switch cols/rows
            [x2 y2] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(end));
        else
            x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
            y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
        end
        %       x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
        %       y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
    case 2
        x1=Args.NumericArguments{1};x2=x1;
        y1=Args.NumericArguments{2};y2=y1;
    case 4
        x1=Args.NumericArguments{1};x2=x1+Args.NumericArguments{3}-1;
        y1=Args.NumericArguments{2};y2=y1+Args.NumericArguments{4}-1;
    otherwise
        error('subaxis argument error')
end


cellwidth=((1-Args.MarginLeft-Args.MarginRight)-(Args.cols-1)*Args.SpacingHorizontal)/Args.cols;
cellheight=((1-Args.MarginTop-Args.MarginBottom)-(Args.rows-1)*Args.SpacingVertical)/Args.rows;
xpos1=Args.MarginLeft+Args.PaddingLeft+cellwidth*(x1-1)+Args.SpacingHorizontal*(x1-1);
xpos2=Args.MarginLeft-Args.PaddingRight+cellwidth*x2+Args.SpacingHorizontal*(x2-1);
ypos1=Args.MarginTop+Args.PaddingTop+cellheight*(y1-1)+Args.SpacingVertical*(y1-1);
ypos2=Args.MarginTop-Args.PaddingBottom+cellheight*y2+Args.SpacingVertical*(y2-1);

if Args.Holdaxis
    h=axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
else
    h=subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
end


set(h,'box','on');
%h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
set(h,'units',get(gcf,'defaultaxesunits'));
set(h,'tag','subaxis');



if (nargout==0), clear h; end;
end
function ArgStruct=parseArgs(args,ArgStruct,varargin)
% Helper function for parsing varargin.
%
%
% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
%
% * ArgStruct is the structure full of named arguments with default values.
% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
% * Aliases can be used to map one argument-name to several argstruct fields
%
%
% example usage:
% --------------
% function parseargtest(varargin)
%
% %define the acceptable named arguments and assign default values
% Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]);
%
% %The capital letters define abrreviations.
% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0)
%
% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
% disp(Args)
%
%
%
%
% Aslak Grinsted 2004

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2double(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (length(varargin)>0)
    FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
    if length(varargin)>1
        Aliases=varargin{2};
    end
end


%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2}); %#ok
FnamesAbbr=strvcat(Fnames{:,3}); %#ok

if length(Aliases)>0
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx)
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not?
        end
        Aliases{i,2}=FieldIdx;
        Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
end

%--------------get parameters--------------------
l=NumArgCount+1;
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx)
            FieldIdx=strmatch(a,FnamesFull);
        end
        if (length(FieldIdx)>1) %shortest fieldname should win
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2};
        end
        
        if isempty(FieldIdx)
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0;
                    end
                else
                    val=true;
                    paramHasValue=0;
                end
            end
            if matlabver>=6
                ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
            else
                ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
            end
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end
end









