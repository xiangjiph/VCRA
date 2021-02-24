function missing = missingTiles(descriptorfolder,directions)
if nargin<2
    directions={'Z'}; % only check z
end
% get list file
myfile = dir(fullfile(descriptorfolder,'list*files'));
fid=fopen(fullfile(descriptorfolder,myfile.name),'r');
inputfiles = textscan(fid,'%s');inputfiles = inputfiles{1};fclose(fid);

folders = unique(cellfun(@fileparts,inputfiles,'UniformOutput',false),'stable');
missing = zeros(length(folders),length(directions));
iter = 0;
for dire = directions{:}
    iter = iter + 1;
    parfor idx = 1:length(folders) 
        matchfile = fullfile(folders{idx},sprintf('match-%s.mat',dire));
        if ~exist(matchfile,'file')
            missing(idx,iter) = 1;
        end
    end
end
