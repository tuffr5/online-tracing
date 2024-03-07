function direction=get_direction(c, plane, prev_plane, cur_idx, center, direction)
    if plane==prev_plane
        return;
    end
    
    % previous center, current center
    cur_xyz=wrap_xyz(c, prev_plane, cur_idx);
    
    prev_xyz=[];
    tr_len=size(center, 1);
    for i=0:tr_len-1
        if i<50 % only a few points needed for local trend
            prev_xyz=[prev_xyz; wrap_xyz([center(tr_len-i,1), center(tr_len-i,2)], center(tr_len-i,3), center(tr_len-i, 4))];
        else
            break;
        end
    end
    
    % flip
    prev_xyz=flip(prev_xyz);
    
    % return direction
    % cannot go back
    if plane==1
        points=[prev_xyz(:,3);cur_xyz(3)];
        direction=interpolate_direction(points', direction);
    end

    if plane==2
        points=[prev_xyz(:,1);cur_xyz(1)];
        direction=interpolate_direction(points', direction);
    end

    if plane==3
        points=[prev_xyz(:,2);cur_xyz(2)];
        direction=interpolate_direction(points', direction);
    end  
end

function direction=interpolate_direction(points, direction)

    if length(points) < 10 % two few points
        return;
    end
        
    % calculate direction
    dist=points;
    dist(end+1)=0;
    dist=dist(2:end);
    
    dist=dist-points;
    dist=dist(1:end-1);
    
    if sum(dist)< -1
        direction='backward';
    else
        direction='forward';
    end
end