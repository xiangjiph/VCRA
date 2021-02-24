function fun_vis_plot_histogram(data, opt)


if isfield(opt,'vis')
    vis = opt.vis;
else
    vis = true;
end
if ~vis
    fig_handle = figure('visible', 'off');
else
    fig_handle = figure('visible', 'on');
end
ax = axes;
histogram(data, opt.plot_opt{:});
set(ax, 'XGrid', 'on', ...
    'YGrid', 'on', ...
    'FontSize', 14, ...
    'Box', true);
% 
% 
if isfield(opt, 'XLabel')
    xlabel(ax, opt.XLabel, 'Interpreter', 'latex');
end
if isfield(opt, 'XScale')
    set(ax, 'XScale', opt.XScale);
end
if isfield(opt, 'YScale')
    set(ax, 'YScale', opt.YScale);
end
if isfield(opt, 'Title')
    title(ax, opt.Title);
end
if isfield(opt, 'YLabel')
    ylabel(ax, opt.YLabel, 'Interpreter', 'latex');
end

if isfield(opt, 'Resolution')
    Resolution = opt.Resolution;
else
    Resolution = 600;
end
if ~isfolder(opt.Output_folder)
    mkdir(opt.Output_folder);
end
if isfield(opt, 'FileName')
    fileType = opt.FileType; 
    if strcmpi(fileType, 'eps')
        print(hfig, '-depsc2', opt.FilePath);
        vers = version();
        if ~strcmp(vers(1:3), '8.4')
            fixPSlinestyle(opt.FilePath);
        end
    elseif strcmpi(fileType, 'pdf')
        print(fig_handle, '-dpdf', opt.FilePath);
    elseif strcmpi(fileType, 'jpg') || strcmpi(fileType, 'jpeg')
        print(fig_handle, '-djpeg', '-opengl', sprintf('-r%d',Resolution), opt.FilePath);
    elseif strcmpi(fileType, 'png') 
        print(fig_handle, '-dpng', '-opengl', sprintf('-r%d',Resolution), opt.FilePath);
    elseif strcmpi(fileType, 'tiff') 
        print(fig_handle, '-dtiff', '-opengl', sprintf('-r%d',Resolution), opt.FilePath);
    elseif strcmpi(fileType, 'emf')
        print(fig_handle, '-dmeta', sprintf('-r%d',Resolution), opt.FilePath); 
    else
        err = MException('', ...
            '=====> ERROR: File type %s is not supported. ', fileType);
        throw(err);
    end
end
end