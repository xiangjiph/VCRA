function transition_str = fun_analysis_percolation_transition_by_bond_removal(graph_str, ...
    bond_occupancy_p_list, num_simulation, parallelQ)

% Search for the transition interval first
bond_rm_f = 1 - bond_occupancy_p_list;
num_search_f = numel(bond_rm_f);
str_template = struct;
[str_template.num_cc, str_template.largest_cc_num_site, str_template.num_site, ...
    str_template.largest_cc_fraction] = deal(nan(num_search_f, 1));
collect_field_name = fieldnames(str_template);

transition_str = struct;
transition_str.num_simulation = num_simulation;
transition_str.bond_occupancy_p = bond_occupancy_p_list;
transition_str.bond_removal_p = bond_rm_f;
transition_str.num_bonds = graph_str.numedges;
transition_str.num_sites = graph_str.numnodes;
record_cell = cell(num_search_f, 1);
if parallelQ
    parfor iter_f = 1 : num_search_f
        tmp_tic = tic;
        record_cell{iter_f} = fun_analysis_percolation_bond_rm_simulation(...
            graph_str, bond_rm_f(iter_f), num_simulation);
        fprintf('Finish simulationg percolation with link removal fraction %f. Elapsed time is %f seconds.\n', ...
            bond_rm_f(iter_f), toc(tmp_tic));
    end
else
    for iter_f = 1 : num_search_f
        tmp_tic = tic;
        record_cell{iter_f} = fun_analysis_percolation_bond_rm_simulation(...
            graph_str, bond_rm_f(iter_f), num_simulation);
        fprintf('Finish simulationg percolation with link removal fraction %f. Elapsed time is %f seconds.\n', ...
            bond_rm_f(iter_f), toc(tmp_tic));
    end   
    
end
%% Organize the data
for iter_f = 1 : num_search_f
    tmp_data = record_cell{iter_f};
    for iter_field = 1 : numel(collect_field_name)
        tmp_field_name = collect_field_name{iter_field};
        transition_str.avg.(tmp_field_name)(iter_f) = mean(tmp_data.(tmp_field_name), 'omitnan');
        transition_str.std.(tmp_field_name)(iter_f) = std(tmp_data.(tmp_field_name), 'omitnan');
    end
end
end