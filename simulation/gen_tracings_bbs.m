addpath(fullfile('skeleton3d'));
addpath(fullfile('skeletonGraph3d'));

tif_dir='simulated_data/*.tif';
swc_dir='simulated_data/swc_gt/';
bb_dir='simulated_data/bbs/';

tifs=dir(tif_dir);

for i=1:numel(tifs)
    [cellNo, chNo]=parse_name(tifs(i).name);
    load(['simulated_data/sim_Cell', cellNo, '_Ch', chNo, 'generated_tissue.mat']);
    
    fid=fopen([bb_dir, 'Cell', cellNo, '_Ch', chNo, '.txt'], 'w');
    for j=1:max(size(colorMatrix))
        thisVol=volumeLabels{1, j};
        thisConnComp=bwconncomp(thisVol, 26);
        if thisConnComp.NumObjects == 1
            % generate tracing from skeleton
            skel = Skeleton3D(thisVol);
            % to tracing
            [swc, seed_xyz]=build_graph_structure(skel);
            if nnz(skel)==length(swc) && ~isempty(swc)
                % generate bbs of 3 planes
                if length(swc) > 30
                    points=swc(1:30, 3:5);
                else
                    points=swc(:, 3:5);
                end
                [xybb, zybb, xzbb]=volume_to_bbs(thisVol, points, seed_xyz, num2str(j));
                fprintf(fid, xybb);
                fprintf(fid, zybb);
                fprintf(fid, xzbb);
                save_to_swc_file(swc, [swc_dir, 'Cell', cellNo, '_Ch', chNo, '/', num2str(j), '.swc']);
            else
                disp(j);
            end
        else
            disp(j);
        end
    end
    fclose(fid);
end


function [xybb, zybb, xzbb]=volume_to_bbs(thisVol, points, seed_xyz, swc_id)
seed_x=seed_xyz(1);
seed_y=seed_xyz(2);
seed_z=seed_xyz(3);

% get mask 
% thisVol -x*y*z
frame_xy=squeeze(thisVol(:,:,seed_z)); % xy
frame_xy=permute(frame_xy, [2,1]);
frame_zy=squeeze(thisVol(seed_x,:,:)); % zy
frame_xz=squeeze(thisVol(:,seed_y,:)); % xz
frame_xz=permute(frame_xz, [2,1]);

% get bb
xybb=get_frame_bb(frame_xy, [seed_y seed_x]);
zybb=get_frame_bb(frame_zy, [seed_y seed_z]);
xzbb=get_frame_bb(frame_xz, [seed_z seed_x]);

[direction_xy,direction_zy,direction_xz]=get_all_directions(points);

xybb=sprintf('%s %s %d %d %d %d %d %d\n', [swc_id '-xy'], direction_xy, xybb(1:4), 1, seed_z);
zybb=sprintf('%s %s %d %d %d %d %d %d\n', [swc_id '-zy'], direction_zy, zybb(1:4), 2, seed_x);
xzbb=sprintf('%s %s %d %d %d %d %d %d\n', [swc_id '-xz'], direction_xz, xzbb(1:4), 3, seed_y);
end

function [direction_xy,direction_zy,direction_xz]=get_all_directions(points)
    dist=points;
    dist(end+1,:)=[0 0 0];
    dist=dist(2:end,:);
    
    dist=dist-points;
    dist=dist(1:end-1,:);
    
    dist_x=sum(dist(:,1));
    dist_y=sum(dist(:,2));
    dist_z=sum(dist(:,3));
    
    if dist_z<0
        direction_xy='backward';
    else
        direction_xy='forward';
    end
    
    if dist_x<0
        direction_zy='backward';
    else
        direction_zy='forward';
    end
    
    if dist_y<0
        direction_xz='backward';
    else
        direction_xz='forward';
    end
end

function bb=get_frame_bb(frame, seed)
CC=bwconncomp(frame);
id = sub2ind(size(frame), seed(1), seed(2));

if CC.NumObjects == 1
    mask=frame;
else
    % get the only object
    for i=1:CC.NumObjects
        if any(CC.PixelIdxList{i} == id)
            mask=zeros(size(frame));
            mask(CC.PixelIdxList{i})=1;
            break;
        end
    end
end
inds=find(mask(:)); 
[y, x]=ind2sub(size(frame),inds(:));
bb=get_bb(min(x):max(x), min(y):max(y));
end

function [cellNo, chNo]=parse_name(name)
str_splited=split(name, '_');
cellNo=str_splited(2);
cellNo=split(cellNo{1}, 'Cell');
cellNo=cellNo(2);
cellNo=cellNo{1};

remaining=str_splited(3);
remaining_splited=split(remaining, '.tif');

chNo=remaining_splited(1);
chNo=split(chNo{1}, 'Ch');
chNo=chNo(2);
chNo=chNo{1};
end