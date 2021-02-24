function [deletenodes,deleteedges] = findloops(A,subs)
%FINDLOOPS Summary of this function goes here
% 
% [OUTPUTARGS] = FINDLOOPS(INPUTARGS) Explain usage here
% 
% Inputs: 
% 
% Outputs: 
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2017/04/05 17:40:20 $	$Revision: 0.1 $
% Copyright: HHMI 2017

nA = size(A,1);
leafnodes = find(sum(A,2)==1);
regnodes = find(sum(A,2)==2);
juncnodes = find(sum(A,2)>2);
criticalnodes = juncnodes;%union(leafnodes,juncnodes);
[~,sortedcriticalnodes] = sort(criticalnodes);
juncvec = zeros(1,nA);
juncvec(juncnodes) = 1;
regvec = zeros(1,nA);
regvec(regnodes) = 1;
leafvec = zeros(1,nA);
leafvec(leafnodes) = 1;
%%
deletenodes = [];
deleteedges = [];
for idx = 1:nA
    if regvec(idx)
        inds=find(A(:,idx));
        % check if there is a connection
        if A(inds(1),inds(2))
            % delete edge
            deletenodes(end+1) = idx;
        end
    elseif juncvec(idx)
        inds=find(A(:,idx));
        inds = inds(juncvec(inds)>0);
        % other junctions
        if length(inds)==2
            % check loop
            if A(inds(1),inds(2)) 
                e1 = sum(abs(subs(inds(1),:)-subs(inds(2),:)));
                e2 = sum(abs(subs(idx,:)-subs(inds(1),:)));
                e3 = sum(abs(subs(idx,:)-subs(inds(2),:)));
                % delete the longest one
                if e1==e2 & e1==e3
                    % perturbe one of the edges
                    subs(idx,:) = subs(idx,:)+.1;
                    e1 = sum(abs(subs(inds(1),:)-subs(inds(2),:)));
                    e2 = sum(abs(subs(idx,:)-subs(inds(1),:)));
                    e3 = sum(abs(subs(idx,:)-subs(inds(2),:)));
                end
                [~,aha] = max([e1 e2 e3]);
                if aha ==1
                    deleteedges(end+1,:) = inds';
                elseif aha==2
                    deleteedges(end+1,:) = [idx inds(1)];
                else
                    deleteedges(end+1,:) = [idx inds(2)];
                end
            end
        else % %% TODO, add a case for triple loops

        end
    else
        continue
    end
end
if ~isempty(deleteedges)
    [deleteedges] = unique([deleteedges;deleteedges(:,[2 1])],'rows');
end
end














