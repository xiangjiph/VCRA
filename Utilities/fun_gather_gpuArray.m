function varargout = fun_gather_gpuArray(input_cell_array)

num_input = numel(input_cell_array);
varargout = cell(num_input, 1);
for output_idx = 1 : num_input
    varargout{output_idx} = gather(input_cell_array{output_idx});
end

end
