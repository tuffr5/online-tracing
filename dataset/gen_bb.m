% based on gt points, determine plane and therefore we can determine
% cur_index, where the bbox is determined by soft segmentation.

base_path='~/Downloads/experiments/csr-dcf/dataset/gt';
swc_dir=dir(fullfile(base_path, '*.swc'));

swc_dir = natsortfiles(swc_dir);

[imgs,h,w] = imread_big('/home/path/to/Downloads/nTracer_sample_img/nTracer_sample_denoised.tif');
%Height*Width*(CZ) --> Height*Width*C*Z !!!important
imgs = reshape(imgs,h,w,3,[]);
% YXCZ -> ZYXC
imgs=permute(imgs,[4 1 2 3]);

outfile=[base_path '_bb.txt'];
fid=fopen(outfile, 'w');

if fid<=0
    error('Fail to open file to write');
end

swc=1;
while swc<=numel(swc_dir)
    % load swc
    swc_points=load_swc_file(fullfile(base_path, swc_dir(swc).name));
    
    swc_id=split(swc_dir(swc).name, '.');
    swc_id=swc_id{1};
    
    % get plane and current index
    if size(swc_points,1)>50
        points=swc_points(1:50,3:5);
    else
        points=swc_points(:,3:5);
    end
    
    [plane, cur_idx, direction, c]=get_plane_idx(points);
    
    % get bb
    [frame, ok]=get_current_frame(imgs, plane, cur_idx);
    if ~ok
        error('can not get current frame. Check the value.');
    end
    
    [ok,bb,fg]=extract_bb(c, frame);
    
    if ~ok
        fprintf(fid, '%s %s - - - - %d %d\n', swc_id, direction, plane, cur_idx);
    else 
        % save to outfile
        fprintf(fid, '%s %s %d %d %d %d %d %d\n', swc_id, direction, bb(1:4), plane, cur_idx);
    end
    swc = swc + 1;
end

fclose(fid);

function [plane, cur_idx, direction, c]=get_plane_idx(points)
    dist=points;
    dist(end+1,:)=[0 0 0];
    dist=dist(2:end,:);
    
    dist=dist-points;
    dist=dist(1:end-1,:);
    
    dist_x=sum(dist(:,1));
    dist_y=sum(dist(:,2));
    dist_z=sum(dist(:,3));
    
    diff=[dist_x dist_y dist_z];
    
    idx=find(abs(diff)==max(abs(diff(:)))); % assumption cross section change faster
    
    switch idx(1)
        case 1
            plane=2;
            cur_idx=points(1,1);
            c=points(1,[3 2]);
        case 2 
            plane=3;
            cur_idx=points(1,2);
            c=points(1,[1 3]);
        case 3
            plane=1;
            cur_idx=points(1,3);
            c=points(1,[1 2]);
        otherwise
            error('No such case');
    end
    
    if diff(idx)<0
        direction='backward';
    else
        direction='forward';
    end
end


function [ok,bb,fg]=extract_bb(c, frame)

template_size=[30 30];
    % extract masked patch: mask out parts outside image
[img, valid_pixels_mask] = get_patch(frame, c, 1.0, template_size);

if max(img(:)) < 20
    ok=false;
    bb=[0,0,0,0];
    fg=[];
    return;
end
    
[h,w,~]=size(img);

% exclude low intensity
img = bsxfun(@times, img, cast(sum(img, 3) > 20, 'like', img));    
% soft segmentation

% to gray
im=rgb2gray(im2double(img));
[I,~]=otsu(im,2);
I=I>1;% exclude background
BW=edge(im,'canny');

% dilation
se0=strel('line',2,0);
se90=strel('line',2,90);
BW=imdilate(BW,[se0 se90]);

BW2=fill_holes(BW);
fg=(I>0)&(BW2>0);


if sum(fg(:))
    % cut out regions outside from image
    fg=fg.*valid_pixels_mask;
    mask=binarize_softmask(fg)>0;
    % get the fg near the center
    L=bwlabel(mask);
    inds=find(L==L(15, 15));

    if ~isempty(inds)
        [y, x]=ind2sub([h w],inds(:));
        y=y+c(2)-0.5*size(mask,1);
        x=x+c(1)-0.5*size(mask,2);
        x(x<1)=0;x(x>size(frame,2))=size(frame,2);
        y(y<1)=0;y(y>size(frame,1))=size(frame,1);
        
        bb=get_bb(min(x):max(x), min(y):max(y));
        ok=true;
    else
        ok=false;
        bb=[0,0,0,0];
    end  
end
end