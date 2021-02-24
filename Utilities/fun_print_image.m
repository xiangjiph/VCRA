function fun_print_image(fig_handle, output_fp)


Resolution = 300;
[out_folder, out_name, file_type] = fileparts(output_fp);
if isempty(file_type)
    file_type = '.png';
    out_name = sprintf('%s%s', out_name, file_type);
    output_fp = fullfile(out_folder, out_name);
end

if ~isfolder(out_folder)
    mkdir(out_folder);
end

switch file_type    
    case '.eps'
        print(fig_handle, '-depsc', '-painters',  output_fp);
    case '.svg'
        print(fig_handle, '-dsvg', output_fp);
    case '.pdf'
        print(fig_handle, '-dpdf', '-bestfit', output_fp);
    case {'.jpg', '.jpeg'}
        print(fig_handle, '-djpeg', '-opengl', sprintf('-r%d',Resolution), output_fp);
    case '.png'
        print(fig_handle, '-dpng', '-opengl', sprintf('-r%d',Resolution), output_fp);
    case {'.tiff', '.tif'}
        print(fig_handle, '-dtiff', '-opengl', sprintf('-r%d',Resolution), output_fp);
    case '.emf'
        print(fig_handle, '-dmeta', sprintf('-r%d',Resolution), output_fp);
    otherwise
        err = MException('', ...
            '=====> ERROR: File type %s is not supported. ', file_type);
        throw(err);
end
fprintf('Finish writing %s\n', output_fp);
end