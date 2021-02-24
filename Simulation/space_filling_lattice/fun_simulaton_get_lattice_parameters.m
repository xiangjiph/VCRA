function lattice_parameter_str = fun_simulaton_get_lattice_parameters(lattice_name)

switch lattice_name
    case '(10, 3)-a'
        a1 = [-1/2; 1/2; -1/2];
        a2 = [1/2; -1/2; -1/2];
        a3 = [1/2; 1/2; -1/2];
        
        b0 = [0; 0; 0];
        b1 = [0; -1/4; 1/4];
        b2 = [-1/4; 1/4; 0];
        b3 = [1/4; 0; -1/4];
        B = [b0, b1, b2, b3];
        bond_length = sqrt(sum((b1 - b0).^2 ));
    case '(10, 3)-b'
        a1 = [sqrt(3); 0; 0];
        a2 = [0; sqrt(3); 0];
        a3 = [sqrt(3)/2; sqrt(3)/2; 3];
        
        b0 = [0; 0; 0];
        b1 = [0; sqrt(3)/2; 1/2];
        b2 = [0; sqrt(3)/2; 3/2];
        b3 = [-sqrt(3)/2; sqrt(3)/2; 2];
        B = [b0, b1, b2, b3];
        bond_length = sqrt(sum((b1 - b0).^2 ));
    case '(10, 3)-c'
        a1 = [sqrt(3); 0; 0];
        a2 = [sqrt(3)/2; 3/2; 0];
        a3 = [0; 0; 9/2];
        
        b0 = [0; 0; 0];
        b1 = [sqrt(3)/2; 0; 1/2];
        b2 = [sqrt(3)/2; 0; 3/2];
        b3 = [sqrt(3)/4; 3/4; 2];
        b4 = [sqrt(3)/4; 3/4; 3];
        b5 = [0; 0; 7/2];
        B = [b0, b1, b2, b3, b4, b5];
        bond_length = sqrt(sum((b1 - b0).^2 ));
    case '(8, 3)-a'
        a1 = [-5/2; sqrt(3)/6; 2 * sqrt(2/3)];
        a2 = [5/2; sqrt(3)/6; 2 * sqrt(2/3)];
        a3 = [0; 4*sqrt(3)/3; sqrt(2/3)];

        b0 = [0; 0; 0];
        b1 = [0; -sqrt(3)/3; -sqrt(2/3)];
        b2 = [1/2; -5*sqrt(3)/6; -sqrt(2/3)];
        b3 = [3/2; -5 * sqrt(3)/6; -sqrt(2/3)];
        b4 = [2; -sqrt(3)/3; -sqrt(2/3)];
        b5 = [2; 0; 0];


        B = [b0, b1, b2, b3, b4, b5];
        bond_length = sqrt(sum((b1 - b0).^2 ));
    case 'cubic'
        a1 = [1; 0; 0];
        a2 = [0; 1; 0];
        a3 = [0; 0; 1];
        
        b0 = [0; 0; 0];
        B = [b0];
        bond_length = 1;
    otherwise 
        error('Unrecognized');
end
translation_mat = [a1, a2, a3];

lattice_parameter_str.lattice_name = lattice_name;
lattice_parameter_str.base_vectors = B;
lattice_parameter_str.translation_vectors = translation_mat;
lattice_parameter_str.bond_length = bond_length;
end