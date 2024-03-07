function [frame, cur_idx, cur_c, ok]=get_next_frame(imgs, bb, prev_plane, plane, direction, ok)
    [d,h,w,c]=size(imgs);
    
    % for error_free
    frame=[];
    
    if prev_plane==plane
        prev_idx=bb(6);
        cur_c=bb(1:2)+0.5*(bb(3:4));
    else
       center=bb(1:2)+0.5*(bb(3:4));% center as pivot
       if prev_plane==1
           if plane==2
               prev_idx=round(center(1));
               cur_c=[bb(6) center(2)];
           end
           
           if plane==3
               prev_idx=round(center(2));
               cur_c=[center(1) bb(6)];
           end
       end
       
       if prev_plane==2
           if plane==1
               prev_idx=round(center(1));
               cur_c=[bb(6) center(2)];
           end
           
           if plane==3
               prev_idx=round(center(2));
               cur_c=[bb(6) center(1)];
           end
       end
       
       if prev_plane==3
           if plane==1
               prev_idx=round(center(2));
               cur_c=[center(1) bb(6)];
           end
           
           if plane==2
               prev_idx=round(center(1));
               cur_c=[center(2) bb(6)];
           end
       end
    end
    
    if plane==1
        if strcmp(direction,'forward')
            cur_idx=prev_idx+1;
            if cur_idx > d
                fprintf('Warning. Exceed depth %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        else
            cur_idx=prev_idx-1;
            if cur_idx < 1
                fprintf('Warning. Illegal index for depth %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        end
        frame=imgs(cur_idx,:,:,:);
    end

    if plane==2
        if strcmp(direction,'forward')
            cur_idx=prev_idx+1;
            if cur_idx > w
                fprintf('Warning. Exceed width %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        else
            cur_idx=prev_idx-1;
            if cur_idx < 1
                fprintf('Warning. Illegal index for width %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        end
        frame=imgs(:,:,cur_idx,:);
        frame=permute(frame, [2,1,3,4]);
    end

    if plane==3
        if strcmp(direction,'forward')
            cur_idx=prev_idx+1;
            if cur_idx > h
                fprintf('Warning. Exceed height %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        else
            cur_idx=prev_idx-1;
            if cur_idx < 1
                fprintf('Warning. Illegal index for height %d.\n',cur_idx-1);
                ok=false;
                return;
            end
        end
        frame=imgs(:,cur_idx,:,:);
    end
    
    frame=squeeze(frame);
    
end