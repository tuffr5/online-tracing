function plane=get_cross_section(imgs, tracker, bb, xyz_plane, debug_file)
    [D,H,W,~]=size(imgs);
    c = tracker.c;
    mask_bb = tracker.mask_bb;
    
    % get each plane using current center as pivot
    plane=bb(5);
    id1=round(c(1)); % corresponding to width
    id2=round(c(2)); % corresponding to height
    id3=bb(6); % correpsonding to plane
    
    [h, w]=size(tracker.mask);
    h=round(1.5*h/2);
    w=round(1.5*w/2);
    
    switch plane
        case 1
            oks=valid_axis([id1, id2, id3], [W, H, D]);   
        case 2
            oks=valid_axis([id1, id2, id3], [D, H, W]); 
        case 3
            oks=valid_axis([id1, id2, id3], [W, D, H]);  
        otherwise
            error('no such plane %d', plane);
    end
    
    % if out of scope, skip cross section change plane by oks
    if sum(oks) < 3
        return;
    end
    
    if plane==1 %xy
        xy_img=imgs(id3,:,:,:);
        zy_img=imgs(:,:,id1,:);
        xz_img=imgs(:,id2,:,:);
        
        [lfwx, rfwx]=valid_bound(id1,w,W);
        [lfhy, rfhy]=valid_bound(id2,h,H);
        [lfwz, rfwz]=valid_bound(id3,w,D);
        [lfhz, rfhz]=valid_bound(id3,h,D);
        
        xy_img=squeeze(xy_img);
        zy_img=squeeze(zy_img);
        xz_img=squeeze(xz_img);

        % convert to the same form of height width
        zy_img=permute(zy_img, [2 1 3]);

        xy_img_cut=xy_img(lfhy:rfhy, lfwx:rfwx, :);
        zy_img_cut=zy_img(lfhy:rfhy, lfwz:rfwz, :);
        xz_img_cut=xz_img(lfhz:rfhz, lfwx:rfwx, :);
    end

    if plane==2 %zy
        xy_img=imgs(id1,:,:,:);
        zy_img=imgs(:,:,id3,:);
        xz_img=imgs(:,id2,:,:);
        
        [lfwz, rfwz]=valid_bound(id1,w,D);
        [lfhy, rfhy]=valid_bound(id2,h,H);
        [lfwx, rfwx]=valid_bound(id3,w,W);
        [lfhx, rfhx]=valid_bound(id3,h,W);
        
        xy_img=squeeze(xy_img);
        zy_img=squeeze(zy_img);
        xz_img=squeeze(xz_img);

        % convert to the same form of height width
        zy_img=permute(zy_img, [2 1 3]);

        xy_img_cut=xy_img(lfhy:rfhy, lfwx:rfwx, :);
        zy_img_cut=zy_img(lfhy:rfhy, lfwz:rfwz, :);
        xz_img_cut=xz_img(lfwz:rfwz, lfhx:rfhx, :);
    end

    if plane==3 %xz
        xy_img=imgs(id2,:,:,:);
        zy_img=imgs(:,:,id1,:);
        xz_img=imgs(:,id3,:,:);
        
        [lfwx, rfwx]=valid_bound(id1,w,W);
        [lfhz, rfhz]=valid_bound(id2,h,D);
        [lfwy, rfwy]=valid_bound(id3,w,H);
        [lfhy, rfhy]=valid_bound(id3,h,H);
        
        xy_img=squeeze(xy_img);
        zy_img=squeeze(zy_img);
        xz_img=squeeze(xz_img);

        % convert to the same form of height width
        zy_img=permute(zy_img, [2 1 3]);

        xy_img_cut=xy_img(lfhy:rfhy, lfwx:rfwx, :);
        zy_img_cut=zy_img(lfwy:rfwy, lfhz:rfhz, :);
        xz_img_cut=xz_img(lfhz:rfhz, lfwx:rfwx, :);
    end
        
    % get fgs in xy/zy/xz plane using hist_fg 
    % to save computation, only calculate around 2xsize(bb)
    % get bboxes over each plane
    [ok_xy,xy_bb,fg_xy]=extract_fg_bb(xy_img_cut, tracker);
    [ok_zy,zy_bb,fg_zy]=extract_fg_bb(zy_img_cut, tracker);
    [ok_xz,xz_bb,fg_xz]=extract_fg_bb(xz_img_cut, tracker);
    
    % determine to use mask_bb or bb, sometimes mask is not accurate
%     if prod(mask_bb(3:4)) > 2*prod(bb(3:4))
%         bb_size=prod(bb(3:4));
%         bb_hw_ratio=bb(3)/bb(4);
%     else
%     [bb_size, bb_hw_ratio]= get_area_hwratio_diff(tracker.mask_bb, 0, 0);
%     end

    bb_size=tracker.maskbb_size_norm;
    bb_hw_ratio=tracker.maskbb_hwratio_norm;
    
    % return most roundish and consistent one
    if ok_xy
        [diff_xy, diff_hw_ratio_xy]= get_area_hwratio_diff(xy_bb, bb_size, bb_hw_ratio);
    else
        diff_xy=100;
        diff_hw_ratio_xy=100; % any number large
    end
    
    if ok_zy
        [diff_zy, diff_hw_ratio_zy]= get_area_hwratio_diff(zy_bb, bb_size, bb_hw_ratio);
    else
        diff_zy=100;
        diff_hw_ratio_zy=100; % any number large
    end
    
    if ok_xz
        [diff_xz, diff_hw_ratio_xz]= get_area_hwratio_diff(xz_bb, bb_size, bb_hw_ratio);
    else
        diff_xz=100;
        diff_hw_ratio_xz=100; % any number large
    end
    
    % add tolerance    
    diff_hw_ratio=[diff_hw_ratio_xy diff_hw_ratio_zy diff_hw_ratio_xz];
    diff=[diff_xy, diff_zy, diff_xz];
    
    if ~isempty(debug_file)
        fprintf(debug_file, '%.2f %.2f %.2f %.2f %.2f %.2f\n', diff, diff_hw_ratio);
    end
    
    norm_diff = diff / norm(diff); % (0,1)
    norm_diff_hw_ratio = diff_hw_ratio / norm(diff_hw_ratio); % (0,1)
    
    norm_exp_diff = exp(-norm_diff.^2 / (2*tracker.maskbb_size_sigma^2));
    norm_exp_diff_hw_ratio = exp(-norm_diff_hw_ratio.^2 / (2*tracker.maskbb_hwratio_sigma^2));
    
    diff_total=0.5*(norm_exp_diff+norm_exp_diff_hw_ratio);

    [~, id]=max(diff_total);

    if id~=plane && diff_total(plane)/diff_total(id)<0.8 ...
        && diff_total(id)>=tracker.cross_plane_degree
        if size(xyz_plane, 1) > 1
            if id~=xyz_plane(end-1,3)
                plane=id;
            end
        end
    end

end

function oks=valid_axis(vals,hs)
oks=[1 1 1];
% h-- high
    for i=1:3
        if vals(i)<1
            oks(i)=0;
            vals(i)=1;
        end
        if vals(i)>hs(i)
            oks(i)=0;
            vals(i)=hs(i);
        end
    end
end

function [bd_s, bd_e]=valid_bound(val,r,t)
% r--range, t -- threshold
% hight threshold
    bd_s=val-r;
    if bd_s < 1
        bd_s=1;
    end
    
    bd_e=val+r;
    if bd_e>t
        bd_e=t;
    end
end