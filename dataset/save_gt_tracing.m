function save_gt_tracing(points, swc_file)
ids_pool=[];% no repetetions
points_pool=[]; % no repetetions
pids_pool=[];

id=1;
pid=-1;

for j=1:size(points,1)
    % find if in ids_pool
    if ~isempty(points_pool)
        C = intersect(points_pool,points(j,:),'rows');
        if ~isempty(C)
            pid=find(ismember(points_pool, C,'rows'));
            pid=pid(1);
        else
            if ~isempty(points_pool) && norm(points_pool(id-1,:)-points(j,:)) > 10
                % two successive points too far, find the nearest in the pool
                dists=pdist2(points_pool, points(j,:));
                pid=find(dists==min(dists(:)));
                pid=pid(1);
            end
        end
    end

    points_pool=[points_pool; points(j,:)];
    ids_pool=[ids_pool; id];
    pids_pool=[pids_pool; pid];
    id=id+1;
    pid=id-1;
end


swc=[];
for i=1:size(points_pool,1)
    % set to swc
    swc = [swc; [ids_pool(i), 0, points_pool(i,1), points_pool(i,2), points_pool(i,3), 0, pids_pool(i)]];
end

save_to_swc_file(swc, swc_file);
end

