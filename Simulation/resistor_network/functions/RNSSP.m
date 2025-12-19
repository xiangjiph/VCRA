classdef RNSSP < handle
    % Resistor network single segment perturbation
    % Use internal conductance matrix - correspond to fix flux boundary
    % condition

    properties
        RN ResistorNetwork
        % Inverse of the weighted laplace matrix of the resistor (without external source)
        InvL
        InvUWL
        %
    end

    %% Methods for initialization
    methods
        function obj = RNSSP(resistor_network)
            obj.RN = resistor_network;
            obj.compute_pseudoinv_L();

        end

        function compute_pseudoinv_L(obj)
            int_c_mat = full(obj.RN.get_internal_conductance_matrix());
            % Compute pseudo inverse - use the properties of the laplace
            % matrix. Seems to be less accurate than directly computing the
            % pseudoinverse (~1% error for InvL * L), takes about 10
            % seconds
            %             mat_size = size(int_c_mat);
            %             assert(mat_size(1) == mat_size(2));
            %             num_node = mat_size(1);
            %             obj.InvL = inv((int_c_mat + 1/num_node)) - 1/num_node;

            % (~ 0.01% error for InvL * L), takes 240 seconds.
            obj.InvL = pinv(int_c_mat);
        end

        function compute_pseudoinv_UWL(obj)
            % Compute the pseudoinver of the unweighted laplacian matrix
            UWL = full(obj.RN.get_unweighted_laplacian_matrix());
            obj.InvUWL = pinv(UWL);
        end

    end

    %% Methods for analysis
    methods
        function [uAu, vAu] = get_ssp_factors_for_node_pair(obj, n_i, n_j, L)
            if nargin < 4
                L = obj.InvL;
            end
            invL_col_i = L(:, n_i);
            invL_col_j = L(:, n_j);

            uAu = invL_col_i(n_i) + invL_col_j(n_j) - ...
                invL_col_i(n_j) - invL_col_j(n_i);

            vAu = invL_col_i(obj.RN.node_pair_idx(:, 1)) + ...
                invL_col_j(obj.RN.node_pair_idx(:, 2)) - ...
                invL_col_i(obj.RN.node_pair_idx(:, 2)) - ...
                invL_col_j(obj.RN.node_pair_idx(:, 1));
        end

        function Fmn_ab = construct_Fmn_ab(obj, invL)
            % [Fmn_ab]_ij: effect of segment j between node a and b on
            % segment i (between node m and n), i.e., each column is the
            % vector for each edge (j)
            % Order of (a, b) and (m, n) matters.
            if nargin < 2
                invL = obj.InvL;
            end
            rn = obj.RN;
            t_tic = tic;
            num_link = numel(rn.map_int_link_label_to_idx);
            Fmn_ab = nan(num_link, num_link);
            for i = 1 : num_link
                ab_idx = rn.map_int_link_label_to_idx(i);
                if ~isnan(ab_idx)
                    node_ab_idx = rn.node_pair_idx(ab_idx, :);
                    assert(all(node_ab_idx > 0));
                    v_a = invL(:, node_ab_idx(1));
                    v_b = invL(:, node_ab_idx(2));
                    fmn_ab = v_a(rn.node_pair_idx(:, 1)) + ...
                        v_b(rn.node_pair_idx(:, 2)) - ...
                        v_a(rn.node_pair_idx(:, 2)) - ...
                        v_b(rn.node_pair_idx(:, 1));
                    Fmn_ab(rn.link_label, i) = rn.link_conductance .* fmn_ab;
                end
            end
            fprintf("Finish constructing Fmn_ab matrix. Elapsed time is %.2f seconds\n", ...
                toc(t_tic));
        end


        function [bin_ind, bin_val] = bin_neighbor_vessel_by_dist(obj, int_link_label, mid_sub, ...
                dist_bin_edge)
            int_link_sub = mid_sub(int_link_label, :);
            dist_to_neighbor = sqrt(sum((mid_sub(obj.RN.link_label, :) - int_link_sub).^2, 2));
            dist_to_neighbor(obj.RN.link_label == int_link_label) = nan;
            [bin_ind, bin_val] = fun_bin_data_to_idx_list_by_edges(dist_to_neighbor , dist_bin_edge);
        end


        function [fsn_in_dist_bin, bin_val] = get_fsn_stat_in_dist_ranges(obj, fsn, int_link_label, ...
                mid_sub, dist_bin_edge)
            if nargin < 5
                dist_bin_edge = 0 : 25 : 500;
            end
            [tmp_bin_ind, bin_val] = obj.bin_neighbor_vessel_by_dist(int_link_label, ...
                mid_sub, dist_bin_edge);

            fsn_in_dist_bin = fun_analysis_get_basic_statistics_in_bins(...
                fsn, tmp_bin_ind);
            fsn_in_dist_bin = fun_analysis_combine_basic_statistics_strs(fsn_in_dist_bin);
        end

        function fss_list = get_fss_factors(obj, invMat, conductance_vec)
            if nargin < 2
                invMat = obj.InvL;
                conductance_vec = obj.RN.link_conductance;
            end
            ind_aa = sub2ind(size(obj.RN.conductance_matrix), obj.RN.node_pair_idx(:, 1), ...
                obj.RN.node_pair_idx(:, 1));
            ind_bb = sub2ind(size(obj.RN.conductance_matrix), obj.RN.node_pair_idx(:, 2), ...
                obj.RN.node_pair_idx(:, 2));
            ind_ab = sub2ind(size(obj.RN.conductance_matrix), obj.RN.node_pair_idx(:, 1), ...
                obj.RN.node_pair_idx(:, 2));
            uAu_list = invMat(ind_aa) + invMat(ind_bb) - 2 * invMat(ind_ab);
            fss_list = uAu_list .* conductance_vec;
        end

    end

    %% Visualization methods
    methods(Static)
        function [fig_hdl, ax_hdl] = plot_loop_geodesic_length_by_vessel_type(loop_geodesic_length, ...
                vessel_tpye, selected_Q)
            if nargin < 3
                selected_Q = true(size(loop_geodesic_length));
            end
            loop_edge = 0.5 : 1 : 25.5;

            int_art_loop = loop_geodesic_length(selected_Q& vessel_tpye==2);
            int_art_loop = int_art_loop(~isnan(int_art_loop));
            art_loop_stat = fun_analysis_get_basic_statistics_with_edges(int_art_loop, loop_edge);

            int_vein_loop = loop_geodesic_length(selected_Q & vessel_tpye==3);
            int_vein_loop = int_vein_loop(~isnan(int_vein_loop));
            vein_loop_stat = fun_analysis_get_basic_statistics_with_edges(int_vein_loop, loop_edge);

            int_cap_loop = loop_geodesic_length(selected_Q & vessel_tpye==1);
            int_cap_loop = int_cap_loop(~isnan(int_cap_loop));
            cap_loop_stat = fun_analysis_get_basic_statistics_with_edges(int_cap_loop, loop_edge);

            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            plot(ax_hdl, art_loop_stat.hist_bin_val, art_loop_stat.hist_pdf, 'LineWidth', 1);
            hold(ax_hdl, 'on');
            plot(ax_hdl, vein_loop_stat.hist_bin_val, vein_loop_stat.hist_pdf, 'LineWidth', 1);
            plot(ax_hdl, cap_loop_stat.hist_bin_val, cap_loop_stat.hist_pdf, 'LineWidth', 1);
            ax_hdl.XLabel.String = 'Number of branch in the shortest loop';
            ax_hdl.YLabel.String = 'PDF';
            ax_hdl.FontSize = 14;
            leg_hdl = legend(ax_hdl, {'Artery', 'Vein', 'Capillary'});
        end

        function [ax_hdl, line_hdl] = plot_subset_conf_int(x, y_array, selection_cell, ax_hdl)
            if nargin < 4
                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
            end

            num_subset = numel(selection_cell);
            line_hdl = [];
            for iter_subset = 1 : num_subset
                tmp_y = y_array(selection_cell{iter_subset}, :);
                tmp_y_stat = fun_analysis_get_basic_statistics_in_column(tmp_y);
                if ~isempty(tmp_y_stat.prctile_val)
                    [ax_hdl, line_hdl(end+1), patch_hdl] = fun_vis_confidence_interval_shaded(...
                        x, tmp_y_stat.median, tmp_y_stat.prctile_val(:, 7), ...
                        tmp_y_stat.prctile_val(:, 9), ax_hdl);
                    patch_hdl.FaceAlpha = 0.2;
                end
                hold(ax_hdl, 'on');
            end
            grid(ax_hdl, 'on');

        end

        function [fig_hdl] = plot_fba_2_1mfaa_vs_euclidean_dist_by_vsl_type(edge_mid_dist, fba_2_1mfaa, ...
                x_bin, y_bin, link_vsl_type, is_internal_link_Q, loglog_fit_min_x)
            fig_hdl = figure;
            ax_hdl_2 = subplot(2,2,4);
            plt_x = edge_mid_dist;
            plt_y = abs(fba_2_1mfaa);
            max_dist = x_bin(end);
            vsl_type = {'Capillary', 'Artery', 'Vein'};
            for i = 1 : 3
                tmp_selected_Q = (link_vsl_type == i) & is_internal_link_Q;
                tmp_x = plt_x(:, tmp_selected_Q);
                tmp_y = plt_y(:, tmp_selected_Q);

                tmp_x_in_range_Q = (tmp_x <= max_dist) & (tmp_x > 0);
                tmp_y = tmp_y(tmp_x_in_range_Q);
                tmp_x = tmp_x(tmp_x_in_range_Q);

                ax_hdl = subplot(2,2,i);
                histogram2(ax_hdl, tmp_x, tmp_y, x_bin, y_bin, ...
                    'DisplayStyle', 'tile');
                hold(ax_hdl, 'on');
                % Overlay in-bin stat
                [bin_idx, bin_val] = fun_bin_data_to_idx_list_by_edges(tmp_x, x_bin);
                bin_stat = fun_analysis_get_basic_statistics_in_bins(tmp_y, bin_idx);
                is_valid_Q = ~arrayfun(@(x) isempty(x.num_data), bin_stat);
                bin_stat_prctile = cat(1, bin_stat(is_valid_Q).prctile_val);
                tmp_eb_x = bin_val(is_valid_Q);
                tmp_eb_y = bin_stat_prctile(:, 8);
                tmp_fit_x_Q = (tmp_eb_x > loglog_fit_min_x);

                errorbar(ax_hdl, tmp_eb_x, tmp_eb_y, ...
                    bin_stat_prctile(:, 8) - bin_stat_prctile(:, 7), ...
                    bin_stat_prctile(:, 9) - bin_stat_prctile(:, 8), 'k', 'LineWidth', 1.5);
                ax_hdl.ColorScale = 'log';
                cbar_hdl = colorbar(ax_hdl);
                cbar_hdl.Limits(1) = 1;
                cbar_hdl.Label.String = 'Number of Data Points';
                tmp_fit_obj = fun_vis_linear_fit_data_in_loglog(tmp_eb_x(tmp_fit_x_Q),...
                    tmp_eb_y(tmp_fit_x_Q), ax_hdl);

                ax_hdl.XLabel.String = 'Distance (\mum)';

                eb_hdl = errorbar(ax_hdl_2, tmp_eb_x, tmp_eb_y, ...
                    bin_stat_prctile(:, 8) - bin_stat_prctile(:, 7), ...
                    bin_stat_prctile(:, 9) - bin_stat_prctile(:, 8), 'LineWidth', 1.5);
                hold(ax_hdl_2, 'on');

                [ax_hdl.YScale, ax_hdl_2.YScale] = deal('log');
                [ax_hdl.XScale, ax_hdl_2.XScale] = deal('log');
                [ax_hdl.YLim(1), ax_hdl_2.YLim(1)] = deal(y_bin(1));
                [ax_hdl.YLabel.Interpreter, ax_hdl_2.YLabel.Interpreter] = deal('latex');
                [ax_hdl.YLabel.String, ax_hdl_2.YLabel.String] = deal('$|f^{\beta}_{\alpha} / (1 - f^{\alpha}_{\alpha})|$');
                ax_hdl.Title.String = sprintf('%s fit min dist %d \\mum', vsl_type{i}, ...
                    loglog_fit_min_x);
            end
            legend(ax_hdl_2, 'Capillary', 'Artery', 'Vein');
            grid(ax_hdl_2, 'on');
            ax_hdl_2.XLabel.String = 'Distance (\mum)';
        end

        function fig_hdl = plot_fba_2_1mfaa_vs_geodesic_dist_by_vsl_type(edge_g_dist, fba_2_1mfaa, ...
                x_bin, y_bin, link_vsl_type, is_internal_link_Q)

            max_g_dist = x_bin(end);
            fig_hdl = figure;
            % fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 2;
            ax_hdl_2 = subplot(2,2,4);
            for i = 1 : 3
                tmp_selected_Q = (link_vsl_type == i) & is_internal_link_Q;
                tmp_x = single(edge_g_dist(:, tmp_selected_Q));
                tmp_y = fba_2_1mfaa(:, tmp_selected_Q);

                tmp_x_in_range_Q = (tmp_x <= max_g_dist) & (tmp_x > 0) & ...
                    isfinite(tmp_x) & isfinite(tmp_y);
                tmp_y = tmp_y(tmp_x_in_range_Q);
                tmp_x = tmp_x(tmp_x_in_range_Q);

                ax_hdl = subplot(2,2,i);
                tmp_hist_hdl = histogram2(ax_hdl, tmp_x, tmp_y, x_bin, y_bin, ...
                    'DisplayStyle', 'tile');
                hold(ax_hdl, 'on');
                % Overlay in-bin stat
                [bin_idx, bin_val] = fun_bin_data_to_idx_list_by_edges(tmp_x, x_bin);
                bin_stat = fun_analysis_get_basic_statistics_in_bins(tmp_y, bin_idx);
                is_valid_Q = ~arrayfun(@(x) isempty(x.num_data), bin_stat);
                bin_stat_prctile = cat(1, bin_stat(is_valid_Q).prctile_val);

                tmp_eb_x = bin_val(is_valid_Q);
                tmp_eb_y = bin_stat_prctile(:, 8);

                errorbar(ax_hdl, tmp_eb_x, tmp_eb_y, ...
                    bin_stat_prctile(:, 8) - bin_stat_prctile(:, 7), ...
                    bin_stat_prctile(:, 9) - bin_stat_prctile(:, 8), 'k', 'LineWidth', 1.5);
                ax_hdl.ColorScale = 'log';
                cbar_hdl = colorbar(ax_hdl);
                cbar_hdl.Limits(1) = 1;
                cbar_hdl.Label.String = 'Number of Data Points';
                ax_hdl.XLabel.String = 'Branch order';

                tmp_fit_obj = fun_vis_linear_fit_data_in_linearlog(tmp_eb_x, ...
                    tmp_eb_y, ax_hdl);

                eb_hdl = errorbar(ax_hdl_2, tmp_eb_x, tmp_eb_y, ...
                    bin_stat_prctile(:, 8) - bin_stat_prctile(:, 7), ...
                    bin_stat_prctile(:, 9) - bin_stat_prctile(:, 8), 'LineWidth', 1.5);
                hold(ax_hdl_2, 'on');

                [ax_hdl.YScale, ax_hdl_2.YScale] = deal('log');
                [ax_hdl.XScale, ax_hdl_2.XScale] = deal('linear');
                % [ax_hdl.YLim(1), ax_hdl_2.YLim(1)] = deal(fnmab_bin(1));
                [ax_hdl.XLim, ax_hdl_2.XLim] = deal(x_bin([1, end]));
                [ax_hdl.YLabel.Interpreter, ax_hdl_2.YLabel.Interpreter] = deal('latex');
                [ax_hdl.YLabel.String, ax_hdl_2.YLabel.String] = deal('$|f^{\beta}_{\alpha} / (1 - f^{\alpha}_{\alpha})|$');

                switch i
                    case 1
                        ax_hdl.Title.String = 'Capillary';
                    case 2
                        ax_hdl.Title.String = 'Artery';
                    case 3
                        ax_hdl.Title.String = 'Vein';
                end
            end
            legend(ax_hdl_2, 'Capillary', 'Artery', 'Vein');
            grid(ax_hdl_2, 'on');
            ax_hdl_2.XLabel.String = 'Branch order';


        end
    end
end