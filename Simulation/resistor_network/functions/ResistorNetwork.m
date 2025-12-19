classdef ResistorNetwork 
    
    properties 
        parameters
        %
        map_node_label_to_idx
        map_node_idx_to_label
        map_node_idx_pair_to_link_label
        num_node
        num_link
        % Internal links 
        link_conductance
        link_label % the link label in the vessel_graph, sorted in ascending order 
        map_int_link_label_to_idx % map the input graph link label to link index in the resistor network
        node_pair_idx
        % Source 
        artery_ep_node_idx
        artery_ep_link_label
        artery_conductance
        artery_pressure
        vein_ep_node_idx
        vein_ep_link_label
        vein_conductance
        vein_pressure
        % Equations
        conductance_matrix
        source
    end
    
    
    methods
        function obj = ResistorNetwork(vessel_graph, valid_node_label, parameters)
            obj.parameters = parameters;
            % Conversion between resistor graph node idx and the vessel
            % graph label
            obj.num_node = numel(valid_node_label);
            obj.map_node_idx_to_label = valid_node_label;
            obj.map_node_label_to_idx = nan(vessel_graph.node.num_cc, 1);
            obj.map_node_label_to_idx(valid_node_label) = 1 : obj.num_node;
            
            % Get the connected links between nodes: 
            [obj.link_conductance, obj.link_label, obj.node_pair_idx] = ...
                obj.get_internal_elements(vessel_graph);
            obj.num_link = numel(obj.link_label);
            % map the input graph link label to link index in the resistor
            % network
            obj.map_int_link_label_to_idx = nan(vessel_graph.link.num_cc, 1);
            obj.map_int_link_label_to_idx(obj.link_label) = 1 : numel(obj.link_label);
            
            obj.map_node_idx_pair_to_link_label = sparse(obj.node_pair_idx(:, 1), ...
                obj.node_pair_idx(:, 2), obj.link_label, obj.num_node, obj.num_node);            
            
            obj.conductance_matrix = obj.get_internal_conductance_matrix();
            
            [obj.artery_ep_node_idx, obj.vein_ep_node_idx, ...
                obj.artery_ep_link_label, obj.vein_ep_link_label, ...
                obj.artery_conductance, obj.vein_conductance, ...
                obj.artery_pressure, obj.vein_pressure] = obj.get_external_elements(vessel_graph);
            [obj.conductance_matrix, obj.source] = obj.add_source();
        end
        
        function [conductance_mat] = get_internal_conductance_matrix(obj)                
            conductance_mat = sparse(obj.node_pair_idx(:, 1), obj.node_pair_idx(:, 2), ...
                obj.link_conductance, obj.num_node, obj.num_node);            
            conductance_mat = conductance_mat + conductance_mat.';
            % Add the diagonal term
            diagonal_term = full(sum(conductance_mat, 2));
            diagonal_term = sparse(1:obj.num_node, 1:obj.num_node, diagonal_term, ...
                obj.num_node, obj.num_node);
            conductance_mat = diagonal_term - conductance_mat;            
        end
        
        function [L] = get_unweighted_laplacian_matrix(obj)
            L = sparse(obj.node_pair_idx(:, 1), obj.node_pair_idx(:, 2), ...
                1, obj.num_node, obj.num_node);            
            L = L + L.';
            % Add the diagonal term
            diagonal_term = full(sum(L, 2));
            diagonal_term = sparse(1:obj.num_node, 1:obj.num_node, diagonal_term, ...
                obj.num_node, obj.num_node);
            L = diagonal_term - L;        
        end
        
        function [link_conductance, link_label, node_pair_idx] = get_internal_elements(...
                obj, vessel_graph)
            % Output: 
            %   link_conductance: 
            %   link_label: the label of the selected link in vessel_graph,
            %   sorted in ascending order
            %   node_pair_idx: the label pair of the nodes connected by the
            %   selected link in vessel_graph. 

            node_label = obj.map_node_idx_to_label;
            % link label in the input vessel graph 
            link_label = cat(1, vessel_graph.node.connected_link_label{node_label});
            link_label = unique(link_label, 'stable');
            
            link_node_pair = vessel_graph.link.connected_node_label(link_label, :);
            % Link should connect to two different node label
            is_valid_pair_Q = ~any(link_node_pair == 0, 2) & ...
                (link_node_pair(:, 1) ~= link_node_pair(:, 2));
            
            link_label = link_label(is_valid_pair_Q);
            link_node_pair = link_node_pair(is_valid_pair_Q, :);
            % Check multi-edge
%             [link_node_pair_u, tmp_ia, tmp_ic] = unique(link_node_pair, ...
%                 'rows', 'stable');
%             if numel(tmp_ia) ~= nnz(is_valid_pair_Q)
%                 warning('Exist multiedge');                
%             end
            
            % Need to figure out why the pruning was not complete. 
            node_pair_idx = obj.map_node_label_to_idx(link_node_pair);
            is_valid_pair_Q = all(~isnan(node_pair_idx), 2);
            node_pair_idx = node_pair_idx(is_valid_pair_Q, :);
            link_label = link_label(is_valid_pair_Q);
%             node_pair_idx = node_pair_idx(is_valid_pair_Q, :);
            % Conductance of multiple links between two nodes are added
            % automatically in the sparse function.
            link_conductance = vessel_graph.link.features.conductance(link_label);
        end
        
        function [artery_ep_node_idx, vein_ep_node_idx, ...
                artery_ep_link_label, vein_ep_link_label, ...
                artery_conductance, vein_conductance, ...
                artery_pressure, vein_pressure] = get_external_elements(obj, vessel_graph)
            % Compute the source part of the conductance matrix and the inhomgenous term
            % This part can be customized later, for different artery and vein pressure
            % Parameters
            artery_label = obj.parameters.artery_label;
            vein_label = obj.parameters.vein_label;
            vein_endpoint_pressure = obj.parameters.vein_endpoint_pressure;
            artery_endpoint_pressure = obj.parameters.artery_endpoint_pressure;
            % Determine the nodes and links directly connect to external source
            link_has_1_endpoint_Q = xor(vessel_graph.link.connected_node_label(:, 1), vessel_graph.link.connected_node_label(:, 2));
            artery_ep_link_label = find(vessel_graph.link.features.type == artery_label & ...
                link_has_1_endpoint_Q);
            vein_ep_link_label = find(vessel_graph.link.features.type == vein_label & ...
                link_has_1_endpoint_Q);
            % Check if the endpoints are in the current connected component(s) or not
            artery_ep_node_label = vessel_graph.link.connected_node_label(artery_ep_link_label, 1);
            vein_ep_node_label = vessel_graph.link.connected_node_label(vein_ep_link_label, 1);
            
            artery_ep_node_idx = obj.map_node_label_to_idx(artery_ep_node_label);
            vein_ep_node_idx = obj.map_node_label_to_idx(vein_ep_node_label);
            
            artery_is_in_cc_Q = ~isnan(artery_ep_node_idx);
            vein_is_in_cc_Q = ~isnan(vein_ep_node_idx);
            
            artery_ep_node_idx = artery_ep_node_idx(artery_is_in_cc_Q);
            vein_ep_node_idx = vein_ep_node_idx(vein_is_in_cc_Q);
            
            artery_ep_link_label = artery_ep_link_label(artery_is_in_cc_Q);
            vein_ep_link_label = vein_ep_link_label(vein_is_in_cc_Q);
            
            % Set the artery and vein endpoint pressure
            if isscalar(artery_endpoint_pressure)
                artery_pressure = ones(size(artery_ep_link_label)) .* artery_endpoint_pressure;
            else
                error('To be implemented');
            end
            if isscalar(vein_endpoint_pressure)
                vein_pressure = ones(size(vein_ep_link_label)) .* vein_endpoint_pressure;
            else
                error('To be implemented');
            end
            artery_conductance = vessel_graph.link.features.conductance(artery_ep_link_label);
            vein_conductance = vessel_graph.link.features.conductance(vein_ep_link_label);
            
        end
        
        function [conductance_mat, source_term] = add_source(obj)
            source_idx = cat(1, obj.artery_ep_node_idx, obj.vein_ep_node_idx);
            source_conductance = cat(1, obj.artery_conductance, obj.vein_conductance);
            source_pressure = cat(1, obj.artery_pressure, obj.vein_pressure);
            
            source_term = sparse(source_idx, 1, source_pressure .* source_conductance, ...
                obj.num_node, 1);
            conductance_to_source_mat = sparse(source_idx, source_idx, source_conductance, ...
                obj.num_node, obj.num_node);
            conductance_mat = obj.conductance_matrix + conductance_to_source_mat;
        end
    end    
end
