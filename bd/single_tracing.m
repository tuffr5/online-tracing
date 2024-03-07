function swc_to_save=single_tracing(bb, direction, imgs, config)
% config
if config.debug
    fid=fopen(config.debug_file, 'w');
else
    fid=[];
end

% files
fileID = fopen(config.txt_file,'w');

% initial bbox, plane, direction
plane=bb(5); % 'xy'

% initialize tracker
switch plane
    case 1 
        frame=imgs(bb(6),:,:,:);
    case 2
        frame=imgs(:,:,bb(6),:);  
        frame=permute(frame, [2,1,3,4]);
    case 3
        frame=imgs(:,bb(6),:,:);
    otherwise
        error("No such plane %d. Please check the value.\n", plane);
end

[d,h,w,c]=size(imgs);

frame=squeeze(frame);
tracker=create_csr_tracker(frame, bb(1:4));

if config.enable_alpha_mask
    tracker.enable_alpha_mask=true;
else
    tracker.enable_alpha_mask=false;
end

tracker.bb6=bb;
tracker.plane=plane;
tracker.prev_plane=plane;
tracker.prev_c=tracker.c;
tracker.direction=direction;

% add to tracker.branch
tracker.branch{1,1}=tracker;

loop_tracker=tracker;
center={};
branch_num=1;
traced=[];
points_pool=[];
swc_to_save=[];
cid=1;

% preset ok
ok=true;

while ~isempty(loop_tracker.branch)
    % organize tracker
    if config.enable_bifurcation
        loop_tracker.branch=organize_branches(loop_tracker.branch, traced);
    end
    
    if isempty(loop_tracker.branch)
        break;
    end
    
    % get the first tracker in tracker.branch
    tracker=loop_tracker.branch{1,1};
    
    xyz=[];
    bb_to_save=[];
    dim_aware=[];
    
    if branch_num>1

        plane=tracker.bb6(5);
        direction=tracker.direction;

        % sanity check for branch bb (due to padding)
        switch plane
            case 1 %xy
                bb=sanity_check_bb(tracker.bb6, w, h);
            case 2 %zy
                bb=sanity_check_bb(tracker.bb6, d, h);
            case 3 %xz
                bb=sanity_check_bb(tracker.bb6, w, d);
            otherwise
                error("No such plane %d. Please check the value.\n", plane);
        end

        % reinitilize tracker.hist
        [frame, ok]=get_current_frame(imgs, plane, bb(6));
        tracker=reinitialize_tracker(frame, bb(1:4), tracker);
        if config.enable_alpha_mask
            tracker.enable_alpha_mask=true;
        else
            tracker.enable_alpha_mask=false;
        end

        tracker.bb6=bb;
        tracker.direction=direction;

        % add the bifurcation point
        xyz=[tracker.main_c tracker.main_bb6(5) tracker.main_bb6(6)];
    end
    
    while ok % detect failure and stop
        bb_to_save=[bb_to_save; bb];
        xyz=[xyz; tracker.c bb(5) bb(6)];
        traced=[traced; bb]; % to aovid duplicates of branch
        points_pool=[points_pool; round(wrap_xyz(tracker.c, bb(5), bb(6)))];
        
        % for swc saving purpose
        if isfield(tracker, 'main_c') %is branch
            tmp_xyz=round(wrap_xyz(tracker.main_c, tracker.main_bb6(5), tracker.main_bb6(6)));
            if ~isempty(points_pool)
                C=intersect(points_pool, tmp_xyz, 'rows');
                if ~isempty(C)
                    pid=find(ismember(points_pool, C, 'rows'));
                    if max(pid)>4 % account for latency of bifurcation
                        pid=pid-4;
                    end
                    if size(pid,1)>1
                        pid=pid(1);
                    end
                end
            end
        else
            pid=cid-1;
        end
        
        swc_to_save=[swc_to_save; [cid, 0, points_pool(end,1), points_pool(end,2), points_pool(end,3), 0, pid]];
        cid=cid+1;
        
        % delete main_bb as tracing goes
        if isfield(tracker, 'main_c') && ~isempty(bb_to_save)
            tracker=rmfield(tracker, {'main_bb6','main_c'});
        end
        

        % visualization and failure detection
        if config.visualize_tracker
            figure(1); if(size(frame,3)<3), colormap gray; end
            imagesc(uint8(frame))
            hold on;
            rectangle('Position',bb(1:4),'LineWidth',1,'EdgeColor','r');
            
            text(15, 25, strcat('plane:', num2str(tracker.prev_plane), ...
                '-index:', num2str(bb(6)), '-', direction), ...
                'Color','r', 'FontSize', 8, 'FontWeight', 'bold');

            hold off;
            truesize;
            drawnow; 
        end
        
        % detect branch in current frame
        if config.enable_bifurcation
            new_tracker=get_bifurcation_by_reconstruction(frame, tracker, bb, config.visualize_bifurcation);
            % add to the end of tracker.branch
            if ~isempty(new_tracker)
                loop_tracker.branch{1,end+1}=new_tracker;
            end
        end
        
        % decide the next frame based on direction and plane
        [frame, cur_idx, cur_c, ok]=get_next_frame(imgs, bb, tracker.prev_plane, plane, direction, ok);
        if ~ok
            break;
        end

        if tracker.prev_plane~=plane
            tracker.c=cur_c;
        end
        % track next frame
        prev_bb=bb;
        tracker.prev_c=tracker.c;
        [tracker, bb, curr_score]=track_csr_tracker(tracker, frame);
        
        if isempty(bb)
            fprintf('Warning. Too low reponse %.2f, stop tracking.\n', curr_score);
            break;
        end
        
        % if bb is nothing but only background stop
        [patch, ~] = get_patch(frame, tracker.c, 1.0, bb(3:4));
        
        dim_aware=[dim_aware max(patch(:)) < tracker.cutoff_intensity];
           
        f = find(diff([false,dim_aware==1,false])~=0);
        [m,~] = max(f(2:2:end)-f(1:2:end-1));
        
        if m>=3
            fprintf('Warning. Too many dim bboxes in a row, stop tracking.\n');
            
            % remove the last m points in xyz, bb_to_save, traced
            xyz=xyz(1:end-m,:);
            bb_to_save=bb_to_save(1:end-m,:);
            traced=traced(1:end-m,:);
            break;
        end
        
        bb=[bb plane cur_idx];

        % determine cross section
        if config.enable_cross_plane
            tracker.prev_plane=plane;
            plane=get_cross_section(imgs, tracker, bb, xyz, fid);

            % determine direction
            direction=get_direction(tracker.c, plane, tracker.prev_plane, cur_idx, xyz, direction);
            tracker.prev_c=tracker.c;
        end
    end
    
    % add to center cell
    if size(xyz, 1) >= 3
        fprintf(fileID,'Branch: %d \n', branch_num);
        
        for i=1:size(bb_to_save,1)
            % save each tracker bb
            if sum(bb_to_save(i, [1 2 6])>0, 'all') == 3
                fprintf(fileID,'%d %d %d %d %d %d\n', round(bb_to_save(i,:)));
            end
        end
        
        center{branch_num}=xyz;
        branch_num=branch_num+1;
    end
    
    % remove the first branch from tracker.branch
    loop_tracker.branch{1,1}=[];
    loop_tracker.branch=loop_tracker.branch(~cellfun('isempty',loop_tracker.branch));
    
    % rest ok
    ok=true;
end


fclose(fileID);

if config.debug
    fclose(fid);
end

% if visualize_tracker
%     figure;
%     center=center(~cellfun('isempty',center));
%     for j=1:size(center, 2)
%         xyz=[];
%         for i=1:size(center{1,j}, 1)
%             cen=center{1,j};
%             xyz=[xyz; wrap_xyz([cen(i,1), cen(i,2)], cen(i,3), cen(i, 4))];
%         end
%         plot3(xyz(:,1), xyz(:,2), xyz(:,3));
%         hold on;
%     end
%     hold off;
% end
end