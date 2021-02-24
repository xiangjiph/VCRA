function output = fun_add_central_difference_terms_in_dir(simu_array, add_dim)


switch add_dim
    case 1
        output = simu_array(1 : (end - 2), 2 : (end - 1), 2 : (end - 1)) + simu_array(3 : (end    ), 2 : (end - 1), 2 : (end - 1));
    case 2
        output = simu_array(2 : (end - 1), 1 : (end - 2), 2 : (end - 1)) + simu_array(2 : (end - 1), 3 : (end    ), 2 : (end - 1));
    case 3
        output = simu_array(2 : (end - 1), 2 : (end - 1), 1 : (end - 2)) + simu_array(2 : (end - 1), 2 : (end - 1), 3 : (end    ));
end
end