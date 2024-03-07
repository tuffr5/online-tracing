function [ok, bb, fg]=extract_fg_bb(img, tracker)
    if max(img(:)) < tracker.cutoff_intensity
        ok=false;
        bb=[0,0,0,0];
        fg=[];
        return;
    end
    
    [h,w,c]=size(img);
    
    % exclude low intensity
    img = bsxfun(@times, img, cast(sum(img, 3) > tracker.cutoff_intensity, 'like', img));    
    % segmentation
    [fg_p, bg_p] = get_location_prior([1, 1, size(img,2), size(img,1)], ...
        tracker.currentScaleFactor*tracker.base_target_sz, [size(img,2), size(img, 1)]);
    [~, fg, ~] = mex_segment(img, tracker.hist_fg, tracker.hist_bg, tracker.nbins, fg_p, bg_p);
    
    if fg
        % cut out regions outside from image
        fg = binarize_softmask(fg);
        
        if tracker.enable_alpha_mask
            if sum(fg(:)) < 1
                fg(round(0.5*size(seg_patch, 2)), round(0.5*size(seg_patch, 1))) = 1;
            end

            % fg = imdilate(fg, strel('disk', 2));
            % fg as trimap
            matting_img=im2double(img);
            % dilation
            fg_trimap=imdilate(fg, strel('disk', 5));
            fg_trimap=0.5*(fg==0 & fg_trimap>0) + fg;
            fg_trimap=fg_trimap.*(fg_trimap > 0 & mean(img, 3) > tracker.cutoff_intensity);
            % matting
            mask=knn_matting_rgb(matting_img, fg_trimap);
            mask=binarize_softmask(mask);

            if sum(mask(:)) > sum(fg(:)) % probably matting is not ok
                fg=mask;
            end
        end
        
        if mask_normal(fg, tracker.target_dummy_area)
            if tracker.mask_diletation_sz > 0
                D = strel(tracker.mask_diletation_type, tracker.mask_diletation_sz);
                fg = imdilate(fg, D);
            end
        else
            pad_h=ceil(h/2-size(tracker.target_dummy_mask, 1)/2);
            pad_w=ceil(w/2-size(tracker.target_dummy_mask, 2)/2);
            if pad_h < 1; pad_h=0; end
            if pad_w < 1; pad_w=0; end
            target_dummy_mask=padarray(tracker.target_dummy_mask, [pad_h pad_w], 0, 'both');
            fg = imresize(target_dummy_mask, size(fg), 'nearest');
        end
        
        [fg, ok, bb]=get_center_fg(logical(fg), tracker);
    end
end