function exit_code = fun_vis_registration_str(registered_str)

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
scatter3(ax_hdl, registered_str.Fixed_sub(:, 2), registered_str.Fixed_sub(:, 1), ...
    registered_str.Fixed_sub(:, 3));
hold(ax_hdl, 'on');
scatter3(ax_hdl, registered_str.Moving_sub(:, 2), registered_str.Moving_sub(:, 1), ...
    registered_str.Moving_sub(:, 3));
legend(ax_hdl, 'Fixed', 'Moving');



end