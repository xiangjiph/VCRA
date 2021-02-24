function fig_handle = fun_vis_two_point_cloud(sub_1, sub_2)


figure;
scatter3(sub_1(:, 1), sub_1(:, 2), sub_1(:, 3));
hold on
scatter3(sub_2(:, 1), sub_2(:, 2), sub_2(:, 3));
legend('1', '2');
end