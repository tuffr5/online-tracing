function new_branch = organize_branches(branch, traced)
% 1. *      *     2.   **
%     *    *          *  *  
%      *  *         *      *
%       **         *        *
% cell_array

if isempty(traced)
    new_branch=branch;
    return;
end

% filter out possible traced points by area
traced_xyz=[];
for i=1:size(traced, 1)
    traced_c=[traced(i,1)+0.5*traced(i,3) traced(i,2)+0.5*traced(i,4)];
    traced_xyz=[traced_xyz; wrap_xyz(traced_c, traced(i,5), traced(i,6))];
end

branch_xyz=[];
for i=1:size(branch,2)
    branch_c=branch{1,i}.c;
    main_bb6=branch{1,i}.main_bb6;
    branch_xyz=[branch_xyz; wrap_xyz(branch_c, main_bb6(5), main_bb6(6))];
    
end

idx=rangesearch(traced_xyz,branch_xyz,10);
for i=1:size(idx,1)
    ids=idx{i};
    if isempty(ids)
        continue;
    else
        for j=1:size(ids,2)
            if i <= size(branch, 2)
                if traced(ids(j),5)~=branch{1,i}.bb6(5) || round(traced(ids(j),6))~=round(branch{1,i}.bb6(6))
                    continue;
                end
                area=rectint(traced(ids(j),1:4), branch{1,i}.bb6(1:4));
                if area/prod(branch{1,i}.bb6(3:4))>0.4 || area/prod(traced(ids(j),3:4))>0.2
                    branch{1,i}=[];
                end
            else
                continue;
            end
            branch=branch(~cellfun('isempty',branch));
        end
    end
end

if size(branch,2) < 1
    new_branch={};
    return;
end

new_branch={};
main_xyz=[];
branch_xyz=[];

for i=1:size(branch,2)
   main_bb6=branch{1,i}.main_bb6;
   main_c=branch{1,i}.main_c;
   branch_c=branch{1,i}.c;
   
   main_xyz=[main_xyz; wrap_xyz(main_c, main_bb6(5), main_bb6(6))];
   branch_xyz=[branch_xyz; wrap_xyz(branch_c, main_bb6(5), main_bb6(6))]; 
end

%range search
idx=rangesearch(branch_xyz,branch_xyz,5);
s=[];
t=[];

for i=1:size(idx,1)
    ids=idx{i};
    if isempty(ids)
        continue;
    else
        s=[s i*ones(1,size(ids,2))];
        t=[t ids];
    end
end

G=graph(s,t);
% remove self_loops
G=rmedge(G,1:numnodes(G),1:numnodes(G));
bins=conncomp(G);

count=1;
for i=1:size(unique(bins),2)
    ids=find(bins==i);
    
    % sanity check, no sudden change allowed
    % assumption: Branch should main position consistence
    % use main_bb as anchors
    now_branch=branch_xyz(ids,:);
    now_main=main_xyz(ids,:);
    diffusion=now_branch-now_main;
    
    
    % get main qudrant
    qudrant=[];
    if sum(diffusion(:,1)>=0) >= sum(diffusion(:,1)<0)% x>=0
        qudrant=[qudrant 1];
    else%x<0
        qudrant=[qudrant -1];
    end
    
    if sum(diffusion(:,2)>=0) >= sum(diffusion(:,2)<0)% y>=0
        qudrant=[qudrant 1];
    else%y<0
        qudrant=[qudrant -1];
    end
    
    if sum(diffusion(:,3)>=0) >= sum(diffusion(:,3)<0)% z>=0
        qudrant=[qudrant 1];
    else%z<0
        qudrant=[qudrant -1];
    end
    
    [id_to_remove, ~]=find(diffusion.*qudrant<0);
    
    if ~isempty(unique(id_to_remove))
        ids(unique(id_to_remove))=[];
    end
    
    if isempty(ids) || length(ids) < branch{1,ids(1)}.branch_counter
        continue;
    end
    
    % calculate direction
    main_direction=branch{1,ids(1)}.direction;
    now_branch=branch_xyz(ids,:);
    now_main=main_xyz(ids,:);
    now_dist=sum(sqrt(power(now_branch-now_main, 2)),2);
    
    dist=now_dist;
    dist(end+1)=0;
    dist=dist(2:end);
    
    dist=dist-now_dist;
    dist=dist(1:end-1);
    
    if sum(dist)>-2 % either far way or closer
        direction=main_direction;
    else
        if strcmp(main_direction, 'forward')
            direction='backward';
        else
            direction='forward';
        end
    end
    
    % add to new_branch and set direction
    new_branch{1,count}=branch{1,ids(1)}; 
    new_branch{1,count}.direction=direction;

    bb=new_branch{1,count}.bb;
    bb6=new_branch{1,count}.bb6;
    main_bb6=new_branch{1,count}.main_bb6;

    padding = new_branch{1,count}.branch_bb_padding;
    
    if prod(bb(3:4)) >= 2*prod(main_bb6(3:4))
        pad_w=round((2*main_bb6(3)-bb(3)));
        pad_h=round((2*main_bb6(4)-bb(4)));
        tmp_bb=bb+[0 0 pad_w pad_h];
        tmp_bb6=bb6+[0 0 pad_w pad_h 0 0];
        if tmp_bb(3) <=2
            tmp_bb(3)=main_bb6(3);
            tmp_bb6(3)=main_bb6(3);
        end
        if tmp_bb(4) <=2
            tmp_bb(4)=main_bb6(4);
            tmp_bb6(4)=main_bb6(4);
        end
        
        new_branch{1,count}.bb=tmp_bb;
        new_branch{1,count}.bb6=tmp_bb6;
    else
        % pad the bb
        new_branch{1,count}.bb=bb+[-padding -padding 2*padding 2*padding];
        new_branch{1,count}.bb6=bb6+[-padding -padding 2*padding 2*padding 0 0];
    end

    count=count+1;
end


