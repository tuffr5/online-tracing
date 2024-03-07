function new_tracker = get_bifurcation_by_edge_detection(img, tracker, bb6, visualize_bifurcation)
    %return a new instance of tracker
    % set plane
    % set direction
    % set bb6
    % set bb
    % set c

    template_size=2*tracker.template_size;
    % extract masked patch: mask out parts outside image
    [seg_patch, valid_pixels_mask] = get_patch(img, tracker.c, 1.0, template_size);
    
    % to gray
    im=rgb2gray(im2double(seg_patch));
    [I,~]=otsu(im,2);
    I=I>1;% exclude background
    BW=edge(im,'canny');
    
    % dilation
    se0=strel('line',2,0);
    se90=strel('line',2,90);
    BW=imdilate(BW,[se0 se90]);
    
    BW2=fill_holes(BW);
    fg=(I>0)&(BW2>0);

    mask=fg.*valid_pixels_mask;
    mask=binarize_softmask(mask)>0;   
    
    % Remove objects containing fewer than 16 pixels
    mask=bwareaopen(mask,4);

    % find connected components
    CC=bwconncomp(mask,8);
    
    nums=CC.NumObjects;
    
    if nums<2
        new_tracker=[];
        return;
    end
    
    masks={};
    bbs=[];
    hists=[];
    areas=[];
    
    for i=1:nums
        if i==1
            mask_cand=bwareafilt(mask,1);
        else
            mask_cand=bwareafilt(mask,i)-bwareafilt(mask,i-1);
        end
        
        masks{i}=mask_cand;
        
        % get_bb
        inds=find(mask_cand);     
        [y, x]=ind2sub([size(mask,1) size(mask,2)],inds(:));
        y=y+tracker.c(2)-0.5*size(mask,1);
        x=x+tracker.c(1)-0.5*size(mask,2);
        bb=get_bb(min(x):max(x), min(y):max(y));
        bbs=[bbs;bb];
        
        % get_hist
        % object rectangle region (to zero-based coordinates)
        obj_reg=[bb(1), bb(2), bb(1)+bb(3), bb(2)+bb(4)] - [1 1 1 1];

        % extract histograms
        hist=mex_extractforeground(img, obj_reg, tracker.nbins, tracker.cutoff_intensity);
        hists=[hists;hist];
        
        % area
        area=rectint(bb, bb6(1:4));
        areas=[areas; area];
        
    end
    
    masks=cat(3, masks{:});
    masks=permute(masks, [3 1 2]);
    
    bbs=reshape(bbs, nums, []);
    areas=reshape(areas, nums, []);
    
    % remove sereve intersected bbs
    idx=[];
    for i=1:nums
        for j=i:nums
            if i==j
                continue;
            end
            area=rectint(bbs(i,:),bbs(j,:));
            if area/min(prod(bbs(i,3:4)), prod(bbs(j,3:4)))>0.4
                %remove the bigger one
                if prod(bbs(i,3:4))>prod(bbs(j,3:4)) 
                    id=i;
                else
                    id=j;
                end
                idx=[idx id];
            end
        end
    end
    
    % remove too dim boxes
    for i=1:nums
        if max(double(seg_patch).*double(repmat(squeeze(masks(i,:,:)), ...
                [1 1 3])), [], 'all') < tracker.cutoff_intensity
            idx=[idx i];
        end
    end
    
    idx=unique(idx);
    if ~isempty(idx)
        masks(idx,:,:)=[];
        bbs(idx,:,:,:,:)=[];
        areas(idx,:)=[];
        hists(idx,:)=[];
    end
    
    nums=size(bbs,1);
    if nums < 2
        new_tracker=[];
        return;
    end
    
    % find the bounding one
    idx=find(areas==max(areas(:)));
    
    % object rectangle region (to zero-based coordinates)
    obj_reg=[bb6(1), bb6(2), bb6(1)+bb6(3), bb6(2)+bb6(4)] - [1 1 1 1];

    % extract histograms
    hist=mex_extractforeground(img, obj_reg, tracker.nbins, tracker.cutoff_intensity);
        
    sims=[];
    for i=1:nums
        sim=cossim(hist', hists(i,:)');
        sims=[sims;sim];
    end
    
    sims=reshape(sims, nums, []);
%     sims=sims/sims(idx(1));
    sims(idx(1))=0;  
    
    id=find(sims==max(sims(:)));
    
    if visualize_bifurcation
        figure(2);
        imagesc(img);
        hold on;

        rectangle('Position',bb6(1:4),'LineWidth',1,'EdgeColor','r');
        rectangle('Position',bbs(id(1),:),'LineWidth',1,'EdgeColor','y');
        rectangle('Position',bbs(idx(1),:),'LineWidth',1,'EdgeColor','b');

        if sims(id(1))>=0.55
            color='r';
        else
            color='g';
        end

        text(15, 25, strcat('max similarity:', num2str(sims(id(1)))), ...
            'Color',color, 'FontSize', 8, 'FontWeight', 'bold');

        hold off;
        truesize;
        drawnow; 
    end
    
    if sims(id(1)) < 0.55% too different
        new_tracker=[];
        return;
    end
    
    bb=bbs(id(1),:);

    new_tracker=tracker;
    new_tracker.c=0.5*[2*bb(1)+bb(3) 2*bb(2)+bb(4)];
    new_tracker.bb=bb;
    new_tracker.bb6=[bb bb6(5) bb6(6)];
    
    new_tracker.main_bb6=bb6;
    new_tracker.main_c=tracker.c;
end