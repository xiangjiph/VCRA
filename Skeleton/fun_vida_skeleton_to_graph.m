function graph_str = fun_vida_skeleton_to_graph(skeleton_str)
% This function convert the skeleton to graph with nodes and links. 
% Input:
%   skeleton_str: structure with fields
%       linear_idx:  N-by-1 double precision array. Index of the skeleton
%       voxel in the mask
%       mask_size: 3-by-1 double precision array. Size of the mask. 
% Output: 
%    graph_str: structure with fields
%       voxel_info: N-by-3 double precision array. Each row consists [voxel
%       index in the mask, number of neighboring skeleton voxel, type of
%       the voxel]. For type of the voxel, positive number means
%       this voxel is a linking voxel or an endpoint. Voxels with the same
%       positive number are in the same link. Voxels with the same negative
%       number are in the same node (the number of voxels in a node can be
%       larger than 1). Voxels with type 0 are isolated points that does
%       not connected to anyone else.
%       
%       nodeList: N-by-1 cell array. Inside is N double precision M-by-1
%       array recording the row index of the voxel in the VOXEL_INFO (Not 
%       in the mask array!) in each node.
%       
%       segmentList: N-by-1 cell array. Inside is N double precision M-by-1
%       array recording the row index of the voxel in the VOXEL_INFO in
%       each link. 
% 
% Adapt from VidaCenterlineToNodeSegList2d by Philp Tsai. 
% 1. Add handling of the voxels on the edge.
% 2. Reduce redundancy in the algorithm
%% Initialization
voxel_list = skeleton_str.linear_idx;
num_skeleton_voxel = numel(skeleton_str.linear_idx);
if isfield(skeleton_str, 'mask_size')
    mask_size = skeleton_str.mask_size;
elseif isfield(skeleton_str, 'array_size')
    mask_size = skeleton_str.array_size;
end

mask_size_pad = mask_size + 2;
num_block_voxel = prod(mask_size);
num_block_voxel_padded = prod(mask_size_pad);
% 26 neighbors relative indices position array:
% [tmp1, tmp2, tmp3] = ndgrid(1:3);
tmp1 = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
tmp2 = [1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3];
tmp3 = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3]; 
neighboring_relative_ind = sub2ind(mask_size_pad, tmp1(:), tmp2(:), tmp3(:));
neighboring_relative_ind = neighboring_relative_ind - neighboring_relative_ind(14);
neighboring_relative_ind(14) = [];
%% Initialization
% Keep the voxels on the boundary by padding the array with 1 layer of
% zeros on all the sides. 
[pos_1_list, pos_2_list, pos_3_list] = ind2sub(mask_size, voxel_list);
% This is the coordinate of the skeleton voxels in the reverse lookup
% table. The lookup table takes the coordinates and returns the skeleton
% idx in the pixel list
voxel_idx_padded = sub2ind(mask_size_pad, pos_1_list + 1, pos_2_list + 1, pos_3_list + 1); 
num_neighbor_list = zeros(num_skeleton_voxel,1);
voxel_type = zeros(num_skeleton_voxel,1);
num_neighbor_processed_list = zeros(num_skeleton_voxel, 1);
% Create a sparse matrix for neighboring voxel indices searching.
reverseLUT = sparse(voxel_idx_padded, ones(num_skeleton_voxel, 1), ...
    1:num_skeleton_voxel, num_block_voxel_padded, 1);
%% Classify voxels
for endpoint_idx = 1 : num_skeleton_voxel
    num_neighbor_list(endpoint_idx) = nnz(reverseLUT(voxel_idx_padded(endpoint_idx ,1) + ...
        neighboring_relative_ind));
end
list_endpoint_idx_in_list = find(num_neighbor_list == 1);
num_endpoints = length(list_endpoint_idx_in_list);

list_node_idx_in_list = find(num_neighbor_list >2);
num_nodes = length(list_node_idx_in_list);
%% Convert skeleton into graph
%initialize counters and parameters
segmentCounter = 0;
nodeCounter = 0;
growingSegment = zeros(1,1000);
growInd = 1;
segmentList = cell(1,num_block_voxel);
nodeList = cell(1,num_block_voxel);
% Flag for conditioning
stillPointsLeftFlag = 1;
endpoint_start_idx = 1;
node_start_idx = 1;
%% Algorithm
while stillPointsLeftFlag == 1
    %Find a starting point
    start_idx_in_list = 0;
    next_neighbor_idx_in_mask = 0;
    previous_idx_in_list = 0;
    %First look for any remaining unprocessed endpoints, process all the
    %end points before processing node points. 
    for endpoint_idx = endpoint_start_idx: num_endpoints
        tmp_idx_in_list = list_endpoint_idx_in_list(endpoint_idx);% index for the endpoint in the voxel list
        %Unused point since it hasn't been labeled in the cpList yet.
        if voxel_type(tmp_idx_in_list) == 0 
            % Found an unprocessed endPoint
            endpoint_start_idx = endpoint_idx;
            start_idx_in_list = tmp_idx_in_list;
            segmentCounter = segmentCounter +1;
            growingSegment = zeros(1,1000);
            growInd = 1;
            voxel_type(start_idx_in_list) = segmentCounter;
            break % Found the unprocessed endpoint, no longer need to search for a while. 
        end
    end
    
    % If all the end points have been processed, process the node points
    % that haven't been processed yet. 
    if start_idx_in_list == 0
       %Found no more endpoints, look for any not-fully-processed nodes
       for node_idx = node_start_idx : num_nodes
           tmp_idx_in_list = list_node_idx_in_list(node_idx);
           if num_neighbor_processed_list(tmp_idx_in_list) < num_neighbor_list(tmp_idx_in_list) 
               %found a good node to start from
               nodeID = voxel_type(tmp_idx_in_list) * -1;
               node_start_idx = node_idx;
               if nodeID == 0
                   %this is a new node.
                   start_idx_in_list = tmp_idx_in_list;
                   segmentCounter = segmentCounter + 1;
                   growingSegment = zeros(1,1000);%[startPointPointer];
                   growInd = 1;
                   nodeCounter = nodeCounter + 1;
                   nodeList{nodeCounter} = [];
                   voxel_type(start_idx_in_list) = -1 * nodeCounter;
               else%if nodeID == 0,
                   %this is an established node.
                   start_idx_in_list = tmp_idx_in_list;
                   segmentCounter = segmentCounter + 1;
                   growingSegment = zeros(1,1000);%[startPointPointer];
                   growInd = 1;
               end
               break
           end
       end
    end
    %-----
    if start_idx_in_list == 0%nothing left to process. Terminate the while loop
        stillPointsLeftFlag = 0;
    end
	% Up to here, we have found the voxel list index to start tracking
    %-----
    if stillPointsLeftFlag
        continueTrackingFlag = 1;
        end_of_segment_reached_flag = 0;
        end_of_track_and_segment_reached_flag = 0;
        current_idx_in_list = start_idx_in_list;
        while continueTrackingFlag == 1
            % Add the current point to the growingSegment list
            % The first element of growingSegment is the index of the voxel
            % in the pixelList, which is the same as the index of the voxel
            % in cpList. Note that it is NOT the index of the voxel in the
            % block mask. 
            
            % if continue tracking, add the current list idx to
            % growingSegment
            
            growingSegment(growInd) = current_idx_in_list;
            growInd = growInd+1;
            
            tmp_num_point_neighbor = num_neighbor_list(current_idx_in_list);
            if tmp_num_point_neighbor == 2
                currentPointType = 'segPoint';
            elseif tmp_num_point_neighbor == 1
                currentPointType = 'endPoint';
            elseif tmp_num_point_neighbor > 2
                currentPointType = 'nodePoint';
            elseif tmp_num_point_neighbor  == 0
                currentPointType = 'dotPoint';
            else
                error('Error: The number of neighboring voxels should be non-negative value');
            end
            
            % Look for next voxel to add to growingSegment. Search for
            %  unprocessed voxels
            neighbor_idx_list_in_list = reverseLUT(neighboring_relative_ind + voxel_idx_padded(current_idx_in_list));
            % Find the valid neighboring pixel to start tracking
            % previousPointerFlatList is a 26x1 logical array. Disallow
            % backtracking.
            neighbor_idx_list_in_list = full(neighbor_idx_list_in_list(neighbor_idx_list_in_list > 0 & ...
                neighbor_idx_list_in_list ~= previous_idx_in_list));
          
            % Select the unprocessed ones
            neighbor_idx_list_in_list = neighbor_idx_list_in_list(num_neighbor_list(neighbor_idx_list_in_list) > ...
                num_neighbor_processed_list(neighbor_idx_list_in_list));
            
            
            if isempty(neighbor_idx_list_in_list) 
                end_of_track_and_segment_reached_flag = 1;
            else
                %Check that no connection has already been made
                num_to_process_neighbor = numel(neighbor_idx_list_in_list);
%                 segmentList{segmentCounter} = growingSegment(growingSegment>0);
                for tmp_idx = 1 : num_to_process_neighbor
                    alreadyLinkedFlag = 0;
                    neighbor_idx_in_list = neighbor_idx_list_in_list(tmp_idx);
                    %Check that a connection to this point doesn't already exit.
                    tmp_point_type = voxel_type(neighbor_idx_in_list);
                    
                    if tmp_point_type >= 0 
                        % This neighbor is a segment point or an
                        % end point or an uninitiated point
                        next_neighbor_idx_in_mask = neighbor_idx_in_list;
                        break% break out of neighbor_idx = 1 : num_to_process_neighbor loop
                    else
                        % this neighbor point is a node.
                        % check that connection doesn't already exist
                        % Most nodes have three branches, but some nodes
                        % have more. Thus, nodeList cannot be initialized
                        % without redudancy at the begining. 
                        segments_idx_list_attached_to_neighbor = nodeList{-tmp_point_type};
                        
%                         For each segments that attached to neighbor that
%                         is a node, check if the current voxel has been
%                         added to any of the segments of the neighboring
%                         node yet.
                        for seg_idx = segments_idx_list_attached_to_neighbor
                            current_seg_idx_list_in_list = segmentList{seg_idx};                            
                            num_voxels_in_seg = numel(current_seg_idx_list_in_list);
                            % Test if the current index have been attached
                            % to the segments or not
                            q = find(current_seg_idx_list_in_list == neighbor_idx_in_list);
                            q1 = max(1,q-1); 
                            q2 = min(q+1,num_voxels_in_seg);
                            directlyLinkedTestSegmentPoints = current_seg_idx_list_in_list(q1:q2);
                            foundLink = ismember(current_idx_in_list,directlyLinkedTestSegmentPoints);
                            
                            if foundLink
                                alreadyLinkedFlag = 1;
                                break %Stop looking at the other segments 
                            end
                            
                        end
                        
                    end

                    if alreadyLinkedFlag == 0
                        next_neighbor_idx_in_mask = neighbor_idx_in_list;
                        break%break out of unfilledIter loop
                    end
                                   
                    if tmp_idx == num_to_process_neighbor
                        % Have looked at all viable neighbors - no
                        % legitimate links (otherwise, the loop has been
                        % broken above
                        disp('This should not occur!'); 
                        keyboard
                        end_of_track_and_segment_reached_flag = 1;
                    end
                end
            end


            if next_neighbor_idx_in_mask == 0
                %no valid neighbor found.  terminate track
                end_of_track_and_segment_reached_flag = 1;
            else%if nextNeighborPointer == 0
                %NEXT, MARK THE CONNECTION TO THE VALID NEIGHBOR
                %mark the connection as made on BOTH points.
                num_neighbor_processed_list(current_idx_in_list) = num_neighbor_processed_list(current_idx_in_list) + 1;
                num_neighbor_processed_list(next_neighbor_idx_in_mask) = num_neighbor_processed_list(next_neighbor_idx_in_mask) + 1;                
                %Now take care of tagging the current point

                switch currentPointType
                    case{'endPoint'}
                        %This should only occur for a starting point.
                        %Tag this point with the segmentID
                        voxel_type(current_idx_in_list) = segmentCounter;
                    case{'segPoint'}
                        %Tag this point with the segmentID
                        voxel_type(current_idx_in_list) = segmentCounter;
                    case{'nodePoint'}
                        %if it reached this point in the code, this node should have remaining connections.
                        %(if not, it would have broken out as an end-of-segment)
                        %it must be a starting node, or a continuing node
                        %either way, do not increment the nodeCount - that was done either at the start point.
                        %or at the last stage, when this node was the next-neighbor.
                        %so, here, simply attach the current segment to the nodeList
                        currentNodeID = - voxel_type(current_idx_in_list);
                        nodeList{currentNodeID} = [nodeList{currentNodeID},segmentCounter];
                end

                % Next, take care of the next neighbor point
                % First determine what type of point it is
                if num_neighbor_list(next_neighbor_idx_in_mask) == 0, nextNeighborType = 'dotPoint';end
                if num_neighbor_list(next_neighbor_idx_in_mask) == 1, nextNeighborType = 'endPoint';end
                if num_neighbor_list(next_neighbor_idx_in_mask) == 2, nextNeighborType = 'segPoint';end
                if num_neighbor_list(next_neighbor_idx_in_mask) > 2, nextNeighborType = 'nodePoint';end
                
                switch nextNeighborType
                    case{'endPoint'}
                        %The track ends here.  Segment termination dealt with in track-ending code below.
                        %Tag the nextNeighbor point with the segmentID
                        voxel_type(next_neighbor_idx_in_mask) = segmentCounter;
                        end_of_track_and_segment_reached_flag = 1;
                    case{'segPoint'}
                        %Continue onward to next cycle
                        previous_idx_in_list = current_idx_in_list;
                        current_idx_in_list = next_neighbor_idx_in_mask;
                        next_neighbor_idx_in_mask = 0;
                    case{'nodePoint'}
                        %Determine the remaining connections on that node
                        %(this is after the current connection is already ascribed)
                        num_remaining_connections_on_neighbor_node = num_neighbor_list(next_neighbor_idx_in_mask) - ...
                            num_neighbor_processed_list(next_neighbor_idx_in_mask);
                        if num_remaining_connections_on_neighbor_node == 0
                            % This node has been fully processed. Track
                            % ends here. Segment termination dealt with in
                            % track-ending code below
                            end_of_track_and_segment_reached_flag = 1;
                        else
                            % this node has more neighbors. Continue
                            % tracking, but end this segment. Segment
                            % termination dealt with in segment-ending code
                            % below.
                            end_of_segment_reached_flag = 1;
                        end

                        % If the neighbor node is full, add the current segment to its nodeList now 
                        % because the code will break out of the tracking loop after this.
                        nodeID = voxel_type(next_neighbor_idx_in_mask) * -1;
                        if num_remaining_connections_on_neighbor_node == 0
                            %the neighbor node is now full (and must therefore also be established)
                            nodeList{nodeID} = [nodeList{nodeID},segmentCounter];
                        else
                            %the neighbor node is still not full.  Is it new or established?
                            if nodeID == 0
                                %node is new
                                %increment the node count here.
                                %the node count is only incremented 2 places - at a startpoint node
                                %and here when the node is a next-neighbor; not when it is the currentPoint.
                                %initiate this nodeList with the current segment. 
                                nodeCounter = nodeCounter + 1;
                                voxel_type(next_neighbor_idx_in_mask) = -1 * nodeCounter;
                                nodeList{nodeCounter} = segmentCounter;
                            else
                                %node is established
                                %add the current segment to this node
                                nodeList{nodeID} = [nodeList{nodeID},segmentCounter];
                            end
                        end
                end
            end

            
            if end_of_segment_reached_flag == 1
                %Our next step is to a node (end of segment) but we can
                %still track onwards, add the neighbor node to the
                %growingSegment, and send it to the segmentList. increment
                %the segmentCounter initiate a new growingSegment with the
                %neighborPoint
                growingSegment(growInd) = next_neighbor_idx_in_mask;
                segmentList{segmentCounter} = growingSegment(growingSegment>0);
                segmentCounter = segmentCounter+1;
                growingSegment = zeros(1,1000);
                growInd = 1;
                previous_idx_in_list = current_idx_in_list;
                current_idx_in_list = next_neighbor_idx_in_mask;
                next_neighbor_idx_in_mask = 0;
                end_of_segment_reached_flag = 0;
            end


            if end_of_track_and_segment_reached_flag==1
                % Our next step is to an termination of the current track.
                % Add the neighbor node to the growingSegment, and send it
                % to the segmentList. Do NOT increment the segmentCounter
                % (the search for a new start point will increment it)
                % initiate a new growingSegment as empty.(the search for a
                % new start point with seed it)
                growingSegment(growInd) = next_neighbor_idx_in_mask;
                segmentList{segmentCounter} = growingSegment(growingSegment>0);
                growingSegment = zeros(1,1000);
                growInd = 1;
                end_of_track_and_segment_reached_flag = 0;
                continueTrackingFlag = 0;
            end
        end
    end
end
%-------
nodeList = nodeList(1:nodeCounter);
segmentList = segmentList(1:segmentCounter);
graph_str.cpList = cat(2, voxel_list, num_neighbor_list, voxel_type);
% graph_str.radiusList = radiusList;
graph_str.segmentList = segmentList;
graph_str.nodeList = nodeList;
graph_str.maskSize = mask_size;

end
