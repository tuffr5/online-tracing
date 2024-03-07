function new_tracker = get_bifurcation_by_reconstruction(img, tracker, bb6, visualize_bifurcation)
    %return a new instance of tracker
    % set plane
    % set direction
    % set bb6
    % set bb
    % set c

    template_size=2*tracker.template_size;
    % extract masked patch: mask out parts outside image
    [seg_patch, valid_pixels_mask] = get_patch(img, tracker.c, 1.0, template_size);
    
    bb_patch=round(0.5*[size(seg_patch,1)-tracker.bb(3),...
                  size(seg_patch,2)-tracker.bb(4), 2*tracker.bb(3), 2*tracker.bb(4)]);
    
    % segmentation
    [fg_p, bg_p] = get_location_prior([1, 1, size(seg_patch,2), size(seg_patch,1)], ...
        tracker.currentScaleFactor*tracker.base_target_sz, [size(seg_patch,2), size(seg_patch, 1)]);
    
    [~, fg, ~] = mex_segment(seg_patch, tracker.hist_fg, tracker.hist_bg, tracker.nbins, fg_p, bg_p);
    
    fg_prob = fg; % save for Astar purpose
    fg = binarize_softmask(fg);
    if tracker.enable_alpha_mask
        if sum(fg(:)) < 1
            fg(round(0.5*size(seg_patch, 2)), round(0.5*size(seg_patch, 1))) = 1;
        end

        % fg as trimap
        matting_img=normalize_img(im2double(seg_patch));
        % dilation
        fg_trimap=0.5*(mean(seg_patch, 3) > tracker.cutoff_intensity & fg < 1) + fg;
        % matting
        alpha=knn_matting_rgb(matting_img, fg_trimap);
        alpha=alpha.*(alpha > 0 & mean(seg_patch, 3) > tracker.cutoff_intensity);
        alpha=binarize_softmask(alpha);

        mask=alpha.*valid_pixels_mask;
        mask=binarize_softmask(mask)>0;
    else
        mask=fg;
    end   
    
    % Remove objects containing fewer than 4 pixels
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
    means_of_masks=[];
    
    for i=1:nums
        mask_cand=zeros(size(mask));
        mask_cand(CC.PixelIdxList{i})=1;
        
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
        
        % mean
        intensity_of_mask=im2double(seg_patch).*double(mask_cand);
        me=sum(intensity_of_mask, [1 2])/sum(mask_cand(:));
        means_of_masks=[means_of_masks; squeeze(me)'];  
    end
    
    masks=cat(3, masks{:});
    masks=permute(masks, [3 1 2]);
    
    bbs=reshape(bbs, nums, []);
    areas=reshape(areas, nums, []);
    min_dists=[];
    
    % find the bounding one
    idx=find(areas==max(areas(:)));
    idx=idx(1);
    
    to_remove_ids=[];
    % remove too far objects
    mask_idx=squeeze(masks(idx,:,:));
    inds_idx=find(mask_idx);     
    [y_idx, x_idx]=ind2sub([size(mask_idx,1) size(mask_idx,2)],inds_idx(:));
    point_idx=[x_idx, y_idx];
    
    for i=1:nums
        if i == idx
            min_dist=inf;
        else
            mask_tmp=squeeze(masks(i,:,:));
            inds_tmp=find(mask_tmp);     
            [y_tmp, x_tmp]=ind2sub([size(mask_idx,1) size(mask_idx,2)],inds_tmp(:));
            point_tmp=[x_tmp, y_tmp];
            dist_tmp=pdist2(point_idx, point_tmp);
            min_dist=min(dist_tmp(:));
        end
        min_dists=[min_dists, min_dist];
        if min_dist > 5
            to_remove_ids=[to_remove_ids i];
        end
    end
        
    
    % remove sereve intersected bbs 
    for i=1:nums
        for j=i:nums
            if i==j
                continue;
            end
            area=rectint(bbs(i,:),bbs(j,:));
            if area/min(prod(bbs(i,3:4)), prod(bbs(j,3:4)))>0.4
                %remove the bigger one
                if i~=idx
                    id=i;
                elseif j~=idx
                    id=j;
                elseif prod(bbs(i,3:4))>prod(bbs(j,3:4))
                    id=i;  
                else
                    id=j;
                end
                to_remove_ids=[to_remove_ids id];
            end
        end
    end
    
    % remove too dim boxes
    for i=1:nums
        if means_of_masks(i) < tracker.cutoff_intensity/255
            to_remove_ids=[to_remove_ids i];
        end
    end
    
    % remove thin lines and too big mask
    for i=1:nums
        fg=squeeze(masks(i,:,:));
        inds=find(fg);     
        [y, x]=ind2sub([size(fg,1) size(fg,2)],inds(:));
        if max(y)-min(y)==0 || max(x)-min(x)==0
            to_remove_ids=[to_remove_ids i];
        else
            try
                bb=get_bb_smallest([x, y]);
                [area, bb_hw_ratio]= get_area_hwratio_diff(bb, 0, 0);
                if bb_hw_ratio > 6 || area > 3*prod(bb6(3:4))
                    to_remove_ids=[to_remove_ids i];
                end
            catch
                to_remove_ids=[to_remove_ids i];
            end
        end
    end
    
    ids=unique(to_remove_ids);
    ids(ids==idx)=[];
    
    nums=size(bbs,1)-length(ids);
    if nums < 2
        new_tracker=[];
        return;
    end
    
    
    % examine every object by backtracking to current traced one
    % filter out by path cost
    probs=[];
    for i=1:size(bbs,1)
        if idx(1) == i % itself
            prob=-inf;
        elseif ~isempty(find(ids-i==0)) %i is in the to_removed mask
            prob=-inf;
        else
            try
                prob=minimum_path_cost_by_prob(seg_patch, bb_patch, fg_prob, ...
                    squeeze(masks(i,:,:)), tracker.branch_color_sigma, visualize_bifurcation);
            catch
                prob=-inf;
            end
        end
        
        probs=[probs prob];
    end
    
    if max(probs(:)) >= tracker.branch_path_degree
        id=find(probs==max(probs(:)));
    else
        id=[];
    end
    
    if isempty(id)
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

function out = normalize_img(img)
    min_val = min(img(:));
    max_val = max(img(:));
    if (max_val - min_val) > 0
        out = (img - min_val)/(max_val - min_val);
    else
        out = zeros(size(img));
    end
end

function prob=minimum_path_cost_by_prob(img, bb1, fg_prob, mask, sigma, visualize)
% gaussian blur
inds=find(mask);     
[y, x]=ind2sub([size(mask,1) size(mask,2)],inds(:));
bb2=get_bb(min(x):max(x), min(y):max(y));

start=get_max_prob_corrd(fg_prob, bb1);
goal=get_max_prob_corrd(fg_prob, bb2);

img=im2double(img);
MAP=permute(img, [2 1 3]);
MAP=rgb2hsv(MAP);

[path, cst]=Astar(MAP, start, goal);

if ~isempty(path)
    prob=exp(-cst.^2 / (2*sigma^2)); % upper bound is 3, lower bound is 0. Theorem 1.
end

if visualize
	figure(2);
	imagesc(img);
	hold on;
    
    if ~isempty(path)
        plot(path(:,1), path(:,2), 'r');

        text(15, 25, strcat('Branch probability:', sprintf('%.2f', prob)), ...
                'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
    else
        text(15, 25, 'No path found', ...
                'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
    end
    rectangle('Position',bb1,'LineWidth',1,'EdgeColor','r');
	rectangle('Position',bb2,'LineWidth',1,'EdgeColor','b');
    
    plot(start(1), start(2), 'x', 'Markersize', 10, 'Color', 'r');
    plot(goal(1), goal(2), 'x', 'Markersize', 10, 'Color', 'b');

	hold off;
	drawnow; 
end

end

function point=get_max_prob_corrd(prob_map, bb)
% find max prob location inside bb
mask=zeros(size(prob_map));
mask(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1)=1;
prob_map=prob_map.*mask;

inds=find(prob_map==max(prob_map(:)));

if length(inds)>1
    mask=zeros(size(prob_map));
    mask(inds)=1;
    % find the biggest one
    CC=bwconncomp(mask,8);
    if CC.NumObjects == 1
        [y, x]=ind2sub([size(mask,1) size(mask,2)],inds(1));
        point=[x, y];
        return;
    else
        numPixels=cellfun(@numel,CC.PixelIdxList);
        [~,idx]=max(numPixels);
        biggest=CC.PixelIdxList{idx};
        [y, x]=ind2sub([size(mask,1) size(mask,2)],biggest(1));
        point=[x, y];
        return;
    end
else
    [y, x]=ind2sub([size(mask,1) size(mask,2)],inds(1));
    point=[x, y];
end
end