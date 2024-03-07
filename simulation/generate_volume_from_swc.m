swc_dir='../dataset/gt/*.swc';
output_dir='simulated_data/material/';

swcs=dir(swc_dir);
swcs = natsortfiles(swcs);

parpool('local',32)
parfor id=1:numel(swcs)
    [retNo, cellNo]=parse_name(swcs(id).name);
    
    outputname=['simulated_data/material/',retNo, '-', cellNo, '.tif'];
%     fprintf("[%3d/%d] Generate %s-%s material...", id, numel(swcs), retNo, cellNo);
    % generate volume from swc
    stack=zeros(512, 512, 136);
    points=load_swc_file(fullfile(swcs(id).folder, swcs(id).name));
    xyz=cutLoc(points(:,3:5));
    idxes=sub2ind(size(stack), xyz(:,2), xyz(:,1), xyz(:,3));
    stack(idxes)=1;
    
    se=strel('sphere', randi([3 6],1));
    stack=imdilate(stack, se);
    
    for i=1:size(stack,3)
        if i==1
            imwrite(squeeze(stack(:,:,i)), outputname);
        else
            imwrite(squeeze(stack(:,:,i)), outputname,'WriteMode','append')
        end
    end
%     fprintf("Finish %s-%s\n", retNo, cellNo);
end

delete(gcp('nocreate'));

function [retNo, cellNo]=parse_name(name)
str_splited=split(name, '-');
retNo=str_splited(1);
retNo=str2num(retNo{1});
retNo=sprintf('%03d', retNo);

remaining=str_splited(2);
remaining_splited=split(remaining, '.swc');

cellNo=remaining_splited(1);
cellNo=str2num(cellNo{1});
cellNo=sprintf('%03d', cellNo);
end

function points=cutLoc(points)
points(points(:,1)<1, 1)=1;
points(points(:,2)<1, 2)=1;
points(points(:,3)<1, 3)=1;
end