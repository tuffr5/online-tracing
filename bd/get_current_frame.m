function [frame, ok] = get_current_frame(imgs, plane, cur_idx)
[d,h,w,~]=size(imgs);
ok=true;
frame=[];

if cur_idx < 1
    fprintf('Warning. Illegal index %d.\n',cur_idx);
    ok=false;
    return;
end
    
if plane==1
    if cur_idx > d
        fprintf('Warning. Exceed depth %d.\n',cur_idx);
        ok=false;
        return;
    end
    frame=imgs(cur_idx,:,:,:);
end

if plane==2
    if cur_idx > w
        fprintf('Warning. Exceed width %d.\n',cur_idx);
        ok=false;
        return;
    end    
    frame=imgs(:,:,cur_idx,:);
    frame=permute(frame, [2,1,3,4]);
end

if plane==3
    if cur_idx > h
        fprintf('Warning. Exceed height %d.\n',cur_idx);
        ok=false;
        return;
    end    
    frame=imgs(:,cur_idx,:,:);
end

frame=squeeze(frame);
end

