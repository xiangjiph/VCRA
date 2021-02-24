function is_simple = fun_skeleton_is_simple_point_removal(neighbor_mat)

% copy neighbors for labeling
num_center_point = size(neighbor_mat,1);
is_simple = true(num_center_point, 1);


% Delete the central point from the neighbor_array and convert it to uint8.
cube = uint8(neighbor_mat(:, [1:13, 15:27]));

label = zeros(num_center_point, 1, 'int8') +2;

% for all points in the neighborhood
for i = 1:26
    
    idx = cube(:,i) == 1 & is_simple;
    
    if any(idx)
        % start recursion with any octant that contains the point i
        switch( i )
            
            case {1,2,4,5,10,11,13}
                cube(idx,:) = octree_labeling(1, label, cube(idx,:) );
            case {3,6,12,14}
                cube(idx,:) = octree_labeling(2, label, cube(idx,:) );
            case {7,8,15,16}
                cube(idx,:) = octree_labeling(3, label, cube(idx,:) );
            case {9,17}
                cube(idx,:) = octree_labeling(4, label, cube(idx,:) );
            case {18,19,21,22}
                cube(idx,:) = octree_labeling(5, label, cube(idx,:) );
            case {20,23}
                cube(idx,:) = octree_labeling(6, label, cube(idx,:) );
            case {24,25}
                cube(idx,:) = octree_labeling(7, label, cube(idx,:) );
            case 26
                cube(idx,:) = octree_labeling(8, label, cube(idx,:) );
        end

        label(idx) = label(idx)+1;
        del_idx = label >= 4;
        
        if any(del_idx)
            is_simple(del_idx) = false;
        end
    end
end
end
%% Sub function
function cube = octree_labeling(octant_label, label, cube)

% check if there are points in the octant with value 1
if( octant_label==1 )
    
    % set points in this octant to current label
    % and recurseive labeling of adjacent octants
    idx = cube(:,1) == 1;
    if any(idx)
        cube(idx,1) = label(idx);
    end
    
    idx = cube(:,2) == 1;
    if any(idx)
        cube(idx,2) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
    end
    
    idx = cube(:,4) == 1;
    if any(idx)
        cube(idx,4) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
    end
    
    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,10) == 1;
    if any(idx)
        cube(idx,10) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
    end
    
    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end
    
end

if( octant_label==2 )
    
    idx = cube(:,2) == 1;
    if any(idx)
        cube(idx,2) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
    end

    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end

    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end

    idx = cube(:,3) == 1;
    if any(idx)
        cube(idx,3) = label(idx);
    end

    idx = cube(:,6) == 1;
    if any(idx)
        cube(idx,6) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,12) == 1;
    if any(idx)
        cube(idx,12) = label(idx);
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

end

if( octant_label==3 )
    
    idx = cube(:,4) == 1;
    if any(idx)
        cube(idx,4) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
    end

    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end

    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end

    idx = cube(:,7) == 1;
    if any(idx)
        cube(idx,7) = label(idx);
    end

    idx = cube(:,8) == 1;
    if any(idx)
        cube(idx,8) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,15) == 1;
    if any(idx)
        cube(idx,15) = label(idx);
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
    
end

if( octant_label==4 )
    
    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
    end

    idx = cube(:,6) == 1;
    if any(idx)
        cube(idx,6) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,8) == 1;
    if any(idx)
        cube(idx,8) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

    idx = cube(:,9) == 1;
    if any(idx)
        cube(idx,9) = label(idx);
    end

    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

end

if( octant_label==5 )
    
    idx = cube(:,10) == 1;
    if any(idx)
        cube(idx,10) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
    end

    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end

    idx = cube(:,18) == 1;
    if any(idx)
        cube(idx,18) = label(idx);
    end

    idx = cube(:,19) == 1;
    if any(idx)
        cube(idx,19) = label(idx);
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end

    idx = cube(:,21) == 1;
    if any(idx)
        cube(idx,21) = label(idx);
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end

    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

end

if( octant_label==6 )
    
    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
    end

    idx = cube(:,12) == 1;
    if any(idx)
        cube(idx,12) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,19) == 1;
    if any(idx)
        cube(idx,19) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
    end


    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,20) == 1;
    if any(idx)
        cube(idx,20) = label(idx);
    end

    idx = cube(:,23) == 1;
    if any(idx)
        cube(idx,23) = label(idx);
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
 
end

if( octant_label==7 )
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = octree_labeling(1,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
    end

    idx = cube(:,15) == 1;
    if any(idx)
        cube(idx,15) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

    idx = cube(:,21) == 1;
    if any(idx)
        cube(idx,21) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
    end

    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end

    idx = cube(:,24) == 1;
    if any(idx)
        cube(idx,24) = label(idx);
    end
    
    idx = cube(:,25) == 1;
    if any(idx)
        cube(idx,25) = label(idx);
        cube(idx,:) = octree_labeling(8,label(idx),cube(idx,:));
    end
end

if( octant_label==8 )
    
    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = octree_labeling(2,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = octree_labeling(3,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = octree_labeling(5,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = octree_labeling(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,23) == 1;
    if any(idx)
        cube(idx,23) = label(idx);
        cube(idx,:) = octree_labeling(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,25) == 1;
    if any(idx)
        cube(idx,25) = label(idx);
        cube(idx,:) = octree_labeling(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,26) == 1;
    if any(idx)
        cube(idx,26) = label(idx);
    end
end
end