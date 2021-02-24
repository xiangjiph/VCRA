function [outputArgs] = writeYML(params, inds, vecfield)
%WRITEYML Summary of this function goes here
% 
% [OUTPUTARGS] = WRITEYML(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/06/09 18:01:41 $	$Revision: 0.1 $
% Copyright: HHMI 2016
pts = vecfield.control;
big = params.big;
fn = params.outfile;
del = sprintf('  ');
% fn = [fn '.control.yml'];
mkdir(fileparts(fn))
fid = fopen(fn, 'w+');
fprintf(fid, '%s\n', ['path: ' params.root]);
% fprintf(fid, '%s\n', ['path: ' tb.root]);
fprintf(fid, 'tiles:\n');

for i = inds
    % make sure path starts with seperator
    if vecfield.path{i}(1)~='/'
        fprintf(fid, '- path: /%s\n', vecfield.path{i});
    else
        fprintf(fid, '- path: %s\n', vecfield.path{i});
    end
    if 0
    tform = eye(5);
    tform(1:3,1:3) = vecfield.afftile(1:3,1:3,i);
    tform(1:3,5) = vecfield.afftile(1:3,4,i);
    else
        tform = vecfield.tform(:,:,i);
    end
    
    fprintf(fid, '%saabb:\n', del);
    fprintf(fid, '%s%sori: [%d, %d, %d]\n', del, del, round(vecfield.origin(i,:)));
    fprintf(fid, '%s%sshape: [%d, %d, %d]\n', del, del, round(vecfield.sz(i,:)));
    
    fprintf(fid, '%sshape:\n', del);
    fprintf(fid, '%s%stype: u16\n', del, del);
    fprintf(fid, '%s%sdims: [%d, %d, %d, %d]\n', del, del, round(params.ymldims));
    
    if big
        fprintf(fid, '%shomography: [%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, tform(1,:));
        fprintf(fid, '%s%s%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, del,  tform(2,:));
        fprintf(fid, '%s%s%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, del,  tform(3,:));
        fprintf(fid, '%s%s%3.1f, %3.1f, %3.1f, %3.1f, %3.1f, ', del, del, tform(4,:));
        fprintf(fid, '%3.1f, %3.1f, %3.1f, %3.1f, %3.1f]\n', tform(5,:));
        
        fprintf(fid, '%sgrid:\n', del);
        fprintf(fid, '%s%sxlims: [%d', del, del, round(vecfield.xlim_cntrl(1)));
        for j = 2:numel(vecfield.xlim_cntrl);
            fprintf(fid, ', %d', round(vecfield.xlim_cntrl(j)));
        end
        fprintf(fid, ']\n');
        fprintf(fid, '%s%sylims: [%d', del, del, round(vecfield.ylim_cntrl(1)));
        for j = 2:numel(vecfield.ylim_cntrl);
            fprintf(fid, ', %d', round(vecfield.ylim_cntrl(j)));
        end
        fprintf(fid, ']\n');
        fprintf(fid, '%s%szlims: [%d', del, del, round(vecfield.zlim_cntrl(1,i)));
        for j = 2:size(vecfield.zlim_cntrl,1);
            fprintf(fid, ', %d', round(vecfield.zlim_cntrl(j,i)));
        end
        fprintf(fid, ']\n');
        %fprintf(fid, '%s%szlims: [%d, %d]\n', del, del, tb.mincrop(i), tb.maxcrop(i));
    else
        fprintf(fid, '%stransform: [%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, tform(1,:));
        fprintf(fid, '%s%s%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, del,  tform(2,:));
        fprintf(fid, '%s%s%17.16f, %17.16f, %17.16f, %3.1f, %17.8f,\n', del, del,  tform(3,:));
        fprintf(fid, '%s%s%3.1f, %3.1f, %3.1f, %3.1f, %3.1f, ', del, del, tform(4,:));
        fprintf(fid, '%3.1f, %3.1f, %3.1f, %3.1f, %3.1f]\n', tform(5,:));
    end
    
    if big
        fprintf(fid, '%s%scoordinates: [%d, %d, %d,\n', del, del, round(pts(1,:,i)));
        for j = 2:size(pts,1)-1
            fprintf(fid, '%s%s%s%d, %d, %d,\n', del, del, del, round(pts(j,:,i)));
        end
        fprintf(fid, '%s%s%s%d, %d, %d]\n', del, del, del, round(pts(end,:,i)));
    end
end
fclose(fid);

end
