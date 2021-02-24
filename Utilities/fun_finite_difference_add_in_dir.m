function x = fun_finite_difference_add_in_dir(x, y, dir)

switch dir
    case -1
        y = y(1 : (end - 2), 2 : (end - 1), 2 : (end - 1));
    case -2
        y = y(2 : (end - 1), 1 : (end - 2), 2 : (end - 1));
    case -3
        y = y(2 : (end - 1), 2 : (end - 1), 1 : (end - 2));
    case 1
        y = y(3 : (end    ), 2 : (end - 1), 2 : (end - 1));
    case 2
        y = y(2 : (end - 1), 3 : (end    ), 2 : (end - 1));
    case 3
        y = y(2 : (end - 1), 2 : (end - 1), 3 : (end    ));     
end
x = x + y;
clearvars y;
end