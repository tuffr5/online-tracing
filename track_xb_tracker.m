function [tracker, region, curr_score] = track_xb_tracker(tracker, img, iteration)
    %% ------------------- TRACKING PHASE -------------------
    % extract features
    f = get_csr_features(img, tracker.c, tracker.currentScaleFactor, ...
        tracker.template_size, tracker.rescale_template_size, ...
        tracker.cos_win, tracker.feature_type, tracker.w2c, tracker.cell_size);

    if ~tracker.use_channel_weights
        response = real(ifft2(sum(fft2(f).*conj(tracker.H), 3)));
    else
        response_chann = real(ifft2(fft2(f).*conj(tracker.H)));
        response = sum(bsxfun(@times, response_chann, reshape(tracker.chann_w, 1, 1, size(response_chann,3))), 3);
    end
    
    % test target presentence for occlusion
    curr_score = [];
    resp_quality = max(response(:));
    
    if isempty(tracker.resp_budg)
        % in first localization frame response score needs to be added as
        % response normalization score
        % therefore 1 is added into the budget
        tracker.resp_norm = resp_quality;
        tracker.resp_budg(end+1) = 1;
    else
        % normalize current response score and add it to the budget
        tracker.resp_budg(end+1) = resp_quality / tracker.resp_norm;
        % if budget is reached, remove first (the oldest) element
        if numel(tracker.resp_budg) > tracker.resp_budg_sz
            tracker.resp_budg(1) = [];
        end
    end
    
    if numel(tracker.resp_budg) == tracker.resp_budg_sz
        response_budget_mean = mean(tracker.resp_budg);
        curr_quality_norm = resp_quality / tracker.resp_norm;
        curr_score = (response_budget_mean - curr_quality_norm) / curr_quality_norm;
        if curr_score > tracker.detect_failure
            region=[];  % target lost
            return;
        end
    end
    
    % find position of the maximum
    [row, col] = ind2sub(size(response),find(response == max(response(:)), 1));

    % calculate detection-based weights
    if tracker.use_channel_weights
        channel_discr = ones(1, size(response_chann, 3));
        for i = 1:size(response_chann, 3)
            norm_response = normalize_img(response_chann(:, :, i));
            local_maxs_sorted = localmax_nonmaxsup2d(squeeze(norm_response(:, :)));

            if local_maxs_sorted(1) == 0, continue; end
            channel_discr(i) = 1 - (local_maxs_sorted(2) / local_maxs_sorted(1));

            % sanity checks
            if channel_discr(i) < 0.5, channel_discr(i) = 0.5; end
        end
    end

    % subpixel accuracy: response map is smaller than image patch -
    % due to HoG histogram (cell_size > 1)
    v_neighbors = response(mod(row + [-1, 0, 1] - 1, size(response,1)) + 1, col);
    h_neighbors = response(row, mod(col + [-1, 0, 1] - 1, size(response,2)) + 1);
    row = row + subpixel_peak(v_neighbors);
    col = col + subpixel_peak(h_neighbors);

    % wrap around 
    if row > size(response,1) / 2
        row = row - size(response,1);
    end
    if col > size(response,2) / 2
        col = col - size(response,2);
    end

    % displacement
    d = tracker.currentScaleFactor * tracker.cell_size * ...
        (1/tracker.rescale_ratio) * [col - 1, row - 1];
    
    % new object center
    c = tracker.c + d;

    % object bounding-box
    region = [c - tracker.currentScaleFactor * tracker.base_target_sz/2, ...
        tracker.currentScaleFactor * tracker.base_target_sz];

    %do a scale space search aswell
    xs = get_scale_subwindow(img, c([2,1]), tracker.base_target_sz([2,1]), ...
        tracker.currentScaleFactor * tracker.scaleSizeFactors, ...
        tracker.scale_window, tracker.scale_model_sz([2,1]), [], tracker.w2c);
    xsf = fft(xs,[],2);
    % scale correlation response
    scale_response = real(ifft(sum(tracker.sf_num .* xsf, 1) ./ (tracker.sf_den + 1e-2) ));
    recovered_scale = ind2sub(size(scale_response),find(scale_response == max(scale_response(:)), 1));
    %set the scale
    currentScaleFactor = tracker.currentScaleFactor * tracker.scaleFactors(recovered_scale);

    % check for min/max scale
    if currentScaleFactor < tracker.min_scale_factor
        currentScaleFactor = tracker.min_scale_factor;
    elseif currentScaleFactor > tracker.max_scale_factor
        currentScaleFactor = tracker.max_scale_factor;
    end
    % new tracker scale
    tracker.currentScaleFactor = currentScaleFactor;
%     disp(currentScaleFactor);

    % put new object location into the tracker structure
    tracker.c = c;
    tracker.bb = region;
    
    %% ------------------- LEARNING PHASE -------------------
    if tracker.use_segmentation
        % convert image in desired colorspace
        if strcmp(tracker.seg_colorspace, 'rgb')
            seg_img = img;
        elseif strcmp(tracker.seg_colorspace, 'hsv')
            seg_img = rgb2hsv(img);
            seg_img = seg_img * 255;
        else
            error('Unknown colorspace parameter');
        end

        % object rectangle region: subtract 1 because C++ indexing starts with zero
        obj_reg = round([region(1), region(2), region(1)+region(3), region(2)+region(4)]) - [1 1 1 1];

        % extract histograms and update them
        hist_fg = mex_extractforeground(seg_img, obj_reg, tracker.nbins, tracker.cutoff_intensity);
        hist_bg = mex_extractbackground(seg_img, obj_reg, tracker.nbins);
        if sum(isnan(hist_fg)) >= 1
            tracker.hist_fg = tracker.hist_fg;
        else
            tracker.hist_fg = (1-tracker.hist_lr)*tracker.hist_fg + tracker.hist_lr*hist_fg;
        end
        
        if sum(isnan(hist_bg)) >= 1
            tracker.hist_bg = tracker.hist_bg;
        else
            tracker.hist_bg = (1-tracker.hist_lr)*tracker.hist_bg + tracker.hist_lr*hist_bg;
        end
        
        % retain for solving color drift
        tracker.hist_fg = (1-tracker.retain_rate)*tracker.hist_fg + tracker.retain_rate*tracker.src_hist_fg;
        tracker.hist_bg = (1-tracker.retain_rate)*tracker.hist_bg + tracker.retain_rate*tracker.src_hist_bg;

        % extract masked patch: mask out parts outside image
        [seg_patch, valid_pixels_mask] = get_patch(seg_img, tracker.c, ...
            tracker.currentScaleFactor, tracker.template_size);

        % segmentation
        [fg_p, bg_p] = get_location_prior([1, 1, size(seg_patch,2), size(seg_patch,1)], ...
            tracker.currentScaleFactor*tracker.base_target_sz, [size(seg_patch,2), size(seg_patch, 1)]);
        [~, fg, ~] = mex_segment(seg_patch, tracker.hist_fg, tracker.hist_bg, tracker.nbins, fg_p, bg_p);
               
        if tracker.enable_alpha_mask
%             save_to_image('analysis/alpha/', strcat(num2str(iteration), 'seg_patch'), seg_patch);
%             save_to_image('analysis/alpha/', strcat(num2str(iteration), 'color_response'), fg);
            % fg as trimap
            fg = binarize_softmask(fg);
            
%             save_to_image('analysis/alpha/', strcat(num2str(iteration), 'color_binarized'), fg);
%             save_to_image('analysis/alpha/', strcat(num2str(iteration), 'color_overlay'), seg_patch, fg);
            if sum(fg(:)) < 1
                fg(round(0.5*size(seg_patch, 2)), round(0.5*size(seg_patch, 1))) = 1;
            end
            matting_img=im2double(seg_patch);
            % dilation
            fg_trimap=imdilate(fg, strel('disk', 3));
            fg_trimap=0.5*(fg==0 & fg_trimap>0) + fg;
            fg_trimap=fg_trimap.*(fg_trimap > 0 & mean(seg_patch, 3) > tracker.cutoff_intensity);
            % matting
            fg=knn_matting_rgb(matting_img, fg_trimap);
            
%             save_to_image('analysis/alpha/', strcat(num2str(iteration), 'alpha_channel'), fg);
        end
        
        fg=binarize_softmask(fg);
        
%         save_to_image('analysis/alpha/', strcat(num2str(iteration), 'alpha_mask'), fg);
%         save_to_image('analysis/alpha/', strcat(num2str(iteration), 'alpha_overlay'), seg_patch, fg);
        % cut out regions outside from image
        mask = single(fg).*single(valid_pixels_mask);
        mask = binarize_softmask(mask);

        % resize to filter size
        mask = imresize(mask, size(tracker.Y), 'nearest');

        % check if mask is too small (probably segmentation is not ok then)
        if mask_normal(mask, tracker.target_dummy_area)
            if tracker.mask_diletation_sz > 0
                D = strel(tracker.mask_diletation_type, tracker.mask_diletation_sz);
                mask = imdilate(mask, D);
            end
        else
            mask = tracker.target_dummy_mask;
        end

    else     
        mask = tracker.target_dummy_mask;
    end
    
    % save for cross_section determination
%     mask=single(bwareafilt(logical(mask), 1));
    [mask, ~, ~]=get_center_fg(logical(mask), tracker);
    inds=find(mask);
    if ~isempty(inds)
        [y, x]=ind2sub(size(tracker.Y),inds(:));
        y=y+tracker.c(2)-0.5*size(mask,1);
        x=x+tracker.c(1)-0.5*size(mask,2);
        mask_bb=get_bb_smallest([double(x), double(y)]);
    end
    
    tracker.mask = mask;
    tracker.mask_bb = mask_bb;

    % mask_bb history for instable reconstrcution of mask
    [bb_size, bb_hw_ratio]=get_area_hwratio_diff(tracker.mask_bb, 0, 0);
    
    tracker.maskbb_size_budg(end+1) = bb_size;
    tracker.maskbb_hwratio_budg(end+1) = bb_hw_ratio;
    % if budget is reached, remove first (the oldest) element
    if numel(tracker.maskbb_size_budg) > tracker.mask_budg_sz
        tracker.maskbb_size_budg(1) = [];
        tracker.maskbb_hwratio_budg(1) = [];
    end
    
    % with a decay kernel
    weight=Quartickernel(1-(1:numel(tracker.maskbb_size_budg))/numel(tracker.maskbb_size_budg));

    tracker.maskbb_size_norm=sum(tracker.maskbb_size_budg.*weight);
    tracker.maskbb_hwratio_norm=sum(tracker.maskbb_hwratio_budg.*weight);

    % extract features from image
    f = get_csr_features(img, tracker.c, tracker.currentScaleFactor, ...
        tracker.template_size, tracker.rescale_template_size, tracker.cos_win, ...
        tracker.feature_type, tracker.w2c, tracker.cell_size);

    % calcualte new filter - using segmentation mask
    H_new = create_xb_filter(f, tracker.Y, single(mask));

    % calculate per-channel feature weights
    if tracker.use_channel_weights
        w_lr = tracker.weight_lr;
        response = real(ifft2(fft2(f).*conj(H_new)));
        chann_w = max(reshape(response, [size(response,1)*size(response,2), size(response,3)]), [], 1) .* channel_discr;
        chann_w = chann_w / sum(chann_w);
        tracker.chann_w = (1-w_lr)*tracker.chann_w + w_lr*chann_w;
        tracker.chann_w = tracker.chann_w / sum(tracker.chann_w);
    end

    % auto-regresive filter update
    lr = tracker.learning_rate;
    tracker.H = (1-lr)*tracker.H + lr*H_new;

    % make a scale search model aswell
    xs = get_scale_subwindow(img, tracker.c([2,1]), tracker.base_target_sz([2,1]), ...
        tracker.currentScaleFactor * tracker.scaleSizeFactors, ...
        tracker.scale_window, tracker.scale_model_sz([2,1]), [], tracker.w2c);
    % fft over the scale dim
    xsf = fft(xs,[],2);
    new_sf_num = bsxfun(@times, tracker.ysf, conj(xsf));
    new_sf_den = sum(xsf .* conj(xsf), 1);
    % auto-regressive scale filters update
    slr = tracker.scale_lr;
    tracker.sf_den = (1 - slr) * tracker.sf_den + slr * new_sf_den;
    tracker.sf_num = (1 - slr) * tracker.sf_num + slr * new_sf_num;

end  % endfunction


function delta = subpixel_peak(p)
	%parabola model (2nd order fit)
	delta = 0.5 * (p(3) - p(1)) / (2 * p(2) - p(3) - p(1));
	if ~isfinite(delta), delta = 0; end
end  % endfunction

function dupl = duplicate_frames(img, img_prev)
    dupl = false;
    I_diff = abs(single(img) - single(img_prev));
    if mean(I_diff(:)) < 0.5
        dupl = true;
    end
end  % endfunction

function [local_max] = localmax_nonmaxsup2d(response)
    BW = imregionalmax(response);
    CC = bwconncomp(BW);

    local_max = [max(response(:)) 0];
    if length(CC.PixelIdxList) > 1
        local_max = zeros(length(CC.PixelIdxList));
        for i = 1:length(CC.PixelIdxList)
            local_max(i) = response(CC.PixelIdxList{i}(1));
        end
        local_max = sort(local_max, 'descend');
    end
end  % endfunction

function out = normalize_img(img)
    min_val = min(img(:));
    max_val = max(img(:));
    if (max_val - min_val) > 0
        out = (img - min_val)/(max_val - min_val);
    else
        out = zeros(size(img));
    end
end  % endfunction

function val = Epanechnikovkernel(x)
    val = 3/4*(ones(size(x,1), size(x,2))-x.^2);
    % rescale so that add up to 1
    val = val / sum(val);
end


function val = Quartickernel(x)
    val = 15/16*(ones(size(x,1), size(x,2))-x.^2).^2;
    % rescale so that add up to 1
    val = val / sum(val);
end

function val = Uniformkernel(x)
    val = ones(size(x,1), size(x,2));
end

function val = Triangularkernel(x)
    val = ones(size(x,1), size(x,2))-x;
    % rescale so that add up to 1
    val = val / sum(val);
end