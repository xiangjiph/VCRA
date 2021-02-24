function [tree1] = appendTree(tree1,tree2)
%APPENDTREE Summary of this function goes here
% 
% [OUTPUTARGS] = APPENDTREE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/21 13:44:01 $	$Revision: 0.1 $
% Copyright: HHMI 2015

% D1 = tree1.D;
% R1 = tree1.R;
% X1 = tree1.X;
% Y1 = tree1.Y;
% Z1 = tree1.Z;
% 
% D2 = tree2.D;
% R2 = tree2.R;
% X2 = tree2.X;
% Y2 = tree2.Y;
% Z2 = tree2.Z;
%%
[y1,x1] = find(tree1.dA);
[y2,x2] = find(tree2.dA);
dims1 = size(tree1.dA);
dims2 = size(tree2.dA);

% correct tree2 indicies
y2 = y2 + dims1(1);
x2 = x2 + dims1(1);

tree1.dA = sparse([y1;y2],[x1;x2],1,dims1(1)+dims2(1),dims1(2)+dims2(2));
tree1.D = [tree1.D;tree2.D];
tree1.R = [tree1.R;tree2.R];
tree1.X = [tree1.X;tree2.X];
tree1.Y = [tree1.Y;tree2.Y];
tree1.Z = [tree1.Z;tree2.Z];

end
