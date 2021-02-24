function tb = LoadBasicYAML(fn)

fid = fopen(fn, 'r');

foundroot = 0;
while ~foundroot
    dat.path = sscanf(fgetl(fid), 'path: %s');
    if ~isempty(dat.path)
        foundroot=1;
    end
end

str = fscanf(fid, '%s');

fclose(fid);

oritok = regexp(str,     'ori:(.\d*.\d*.\d*.)', 'tokens');
shapetok = regexp(str, 'shape:(.\d+.\d+.\d+.)', 'tokens');
dimstok = regexp(str, 'dims:(.\d*.\d*.\d*.\d*.)', 'tokens');
typetok = regexp(str, 'type:(\w{3})', 'tokens');
pathtok = regexp(str, 'path:(\/.*?\d{5})', 'tokens');
tformtok = regexp(str, 'transform:(\[.*?])', 'tokens');
% tformtok = regexp(str, 'homography:(\[.*?])', 'tokens');

Ntiles = numel(oritok);

for i = 1:Ntiles

    dat.tiles{i}.aabb.ori   = eval(oritok{i}{1});
    dat.tiles{i}.aabb.shape = eval(shapetok{i}{1});
    dat.tiles{i}.shape.dims = eval(dimstok{i}{1});
    
    dat.tiles{i}.path       = pathtok{i}{1};
    dat.tiles{i}.shape.type = typetok{i}{1};
    
    dat.tiles{i}.transform  = eval(tformtok{i}{:});
end

tb.tform = zeros(5,5,Ntiles);
tb.origin = zeros(Ntiles, 3);
tb.sz = zeros(Ntiles, 3);
tb.dims = zeros(Ntiles, 4);
tb.bbox = zeros(8,3,Ntiles);
tb.corners = zeros(8,3,Ntiles);
tb.path = cell(Ntiles,1);
tb.inds = zeros(2,3,Ntiles);

for i = 1:numel(dat.tiles)
   
    tb.tform(:,:,i)= reshape(dat.tiles{i}.transform, 5, 5)';
    tb.origin(i,:) = dat.tiles{i}.aabb.ori;
    tb.sz(i,:)     = dat.tiles{i}.aabb.shape;
    tb.dims(i,:)   = dat.tiles{i}.shape.dims;
    tb.bbox(:,:,i) = getBbox(tb.origin(i,:), tb.sz(i,:));
    tb.corners(:,:,i) = projectCorners(tb.tform(:,:,i), tb.dims(i,:));
    tb.path{i} = dat.tiles{i}.path;
    
end
tb.inds(2,:,:) = tb.dims(:,1:3)';
tb.root = dat.path;
tb.N = Ntiles;

