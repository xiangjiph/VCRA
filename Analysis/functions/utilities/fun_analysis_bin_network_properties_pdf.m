function bin_data_str = fun_analysis_bin_network_properties_pdf(data_cell, bin_edge_cell, bin_val)

% Interpolate the PDF and resample it 
num_cell = numel(data_cell);
interpolate_data = zeros(numel(bin_val), num_cell);
is_valid_Q = false(num_cell, 1);
bin_val_step = bin_val(2) - bin_val(1);
for iter_cell = 1 : num_cell
   tmp_data = data_cell{iter_cell};
   tmp_bin_edge = bin_edge_cell{iter_cell};
   if numel(tmp_bin_edge) < 2
       continue;
   end
   % Assume linear edge
   tmp_edge_step = tmp_bin_edge(2) - tmp_bin_edge(1);
   tmp_data = tmp_data .* tmp_edge_step;
   tmp_cdf = [0, cumsum(tmp_data)];
   % Determine the interpolation position
   tmp_bin_val = [bin_val(1) - bin_val_step/2, bin_val + bin_val_step/2];
   tmp_cdf_interpolate = griddedInterpolant(tmp_bin_edge, tmp_cdf, 'linear', 'nearest');
   tmp_cdf_at_bin = tmp_cdf_interpolate(tmp_bin_val); % Linear interpolation
%    tmp_cdf_at_bin = interp1(tmp_bin_edge, tmp_cdf, tmp_bin_val, 'linear');
   tmp_pdf_at_bin = diff(tmp_cdf_at_bin) ./ bin_val_step;
   
   tmp_valid_bin_Q = (bin_val >= (tmp_bin_edge(1) - tmp_edge_step/2)) & ...
       (bin_val <= (tmp_bin_edge(end) + tmp_edge_step/2));   
   
   interpolate_data(tmp_valid_bin_Q, iter_cell) = tmp_pdf_at_bin(tmp_valid_bin_Q);
   is_valid_Q(iter_cell) = true;
end
interpolate_data = interpolate_data(:, is_valid_Q);
interpolate_data = interpolate_data';

bin_data_str.data = interpolate_data;
bin_data_str.bin = bin_val;

for iter_bin = numel(bin_val) : -1 : 1
    bin_stat(iter_bin) = fun_analysis_get_basic_statistics(interpolate_data(:, iter_bin), true);
end
bin_data_str.bin_stat = bin_stat;
end