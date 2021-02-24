function fun_print_image_in_several_formats(fig_handle, output_fp, format_list)

if nargin < 3
    format_list = {'.png', '.eps'};
end
[out_folder, out_name, ~] = fileparts(output_fp);
for iter_format = 1 : numel(format_list)
    tmp_out_name = sprintf('%s%s', out_name, format_list{iter_format});
    output_fp = fullfile(out_folder, tmp_out_name);
    fun_print_image(fig_handle, output_fp);
end

fig_name = fullfile(out_folder, sprintf('%s_data.fig', out_name));
% save(fig_name, 'fig_handle');
savefig(fig_handle, fig_name, 'compact');
fprintf('Finish writing %s\n', fig_name);
% if ~isempty(fig_handle.UserData)
%     xml_fp = sprintf('%s_plot_data.xlm', fullfile(out_folder, out_name));
%     xml_write(xml_fp, fig_handle.UserData);
%     fprintf('Finish writing figure UserData to %s\n', xml_fp);
% end
end