function [fg, ok, bb]=get_center_fg(mask, tracker)
    bbo=[size(mask,2)/2-tracker.bb(3)/2 size(mask,1)/2-tracker.bb(4)/2 tracker.bb(3) tracker.bb(4)];
    
    mask=bwareaopen(mask,4); %some neurons are extremely small
    % find connected components
    CC=bwconncomp(mask,8);
    
    
    if CC.NumObjects < 1
        ok=false;
        bb=[0,0,0,0;
            0,0,0,0];
        fg=[];
        return;
    end

    areas=[];
    masks={};
    bbs=[];
    nums=CC.NumObjects;

    for i=1:nums
        mask_cand=zeros(size(mask));
        mask_cand(CC.PixelIdxList{i})=1;
        masks{i}=mask_cand;

        % get_bb
        inds=find(mask_cand);     
        [y, x]=ind2sub([size(mask,1) size(mask,2)],inds(:));
        bb=get_bb(min(x):max(x), min(y):max(y));
        bbs=[bbs;bb];

        % area
        area=rectint(bb, bbo);
        areas=[areas; area];
    end

    masks=cat(3, masks{:});
    masks=permute(masks, [3 1 2]);
    bbs=reshape(bbs, nums, []);
    areas=reshape(areas, nums, []);
    
    idx=find(areas==max(areas(:)));
    if length(idx)>1
        s=regionprops(CC,'Centroid');
        centroids = cat(1,s.Centroid);
        [~,id]=min(vecnorm(centroids-[bbo(1)+0.5*bbo(3) bbo(2)+0.5*bbo(4)]));
        id=id(1);
    else
        id=1;
    end
%     bb=bbs(idx(id),:);
    ok=true;
    fg=squeeze(masks(idx(id),:,:));
    
    % get smallest bounding box
    inds=find(fg);     
    [y, x]=ind2sub([size(fg,1) size(fg,2)],inds(:));
    if max(y)-min(y)==0 || max(x)-min(x)==0
        ok=false;
        bb=[0,0,0,0;
            0,0,0,0];
        return;
    end
    
    bb=get_bb_smallest([x, y]);
    
end