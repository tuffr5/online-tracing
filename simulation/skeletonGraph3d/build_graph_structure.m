function [swc, seed_xyz] = build_graph_structure(skel, seed_cordinates)
% find path between root to every end points so that the parent_child
% reletionship can be achieved -> convert skel to swc file format
% Date: Apr 1, 2020



w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

skel2=padarray(skel,[1 1 1]);

list_neuron=find(skel2);

% 26-nh of all canal voxels
nh = logical(pk_get_nh(skel2,list_neuron));

% 26-nh indices of all canal voxels
nhi = pk_get_nh_idx(skel2,list_neuron);

% remove center points itself
% nh(:,14) = [];
% nhi(:,14) = [];

sum_nh = sum(logical(nh),2);

for i=1:length(list_neuron)
    cell_child_node{i, 1} = list_neuron(i);
    tmp = nhi(i, find(nh(i,:)));
    tmp(tmp==list_neuron(i)) = [];
    cell_child_node{i,2}= tmp;
end

% configure adj graph

for j=1:size(list_neuron, 1)
    list_child_node = cell_child_node{j, 2}';
    for jj=1:length(list_child_node)
        if j==1 && jj==1
            node_j = j;
            node_jj = find(list_neuron==list_child_node(jj));
        else
            node_j = [node_j j];
            node_jj = [node_jj find(list_neuron==list_child_node(jj))];
        end
    end
end

adj = sparse(node_j, node_jj, ones(size(node_j, 1),1), size(list_neuron, 1), size(list_neuron, 1));

% consrtuct graph
G = graph(adj);
% delete self_loops and loops between two same nodes
H = simplify(G);

% path search between two end points
ep = find(sum_nh==2);

if nargin == 1
    % Determine which dimension sum will use
    seed_idx = 1;
else
    seed = sub2ind([w+2,l+2,h+2], seed_cordinates(1)+1, seed_cordinates(2)+1, seed_cordinates(3)+1);
    seed_idx = find(list_neuron(ep)==seed);
end

ep_diff = setdiff(1:length(ep), seed_idx);

% from seed_point to other ep
if ~isempty(ep_diff)

    [TR,~] = shortestpathtree(H,ep(seed_idx),ep(ep_diff), 'OutputForm','cell');
    [~,TR_sort_index] = sort(cellfun(@length,TR), 'descend');
    count = 0;
    for i=1:numel(TR)
        % maybe a path per swc
        path = TR{TR_sort_index(i)};
        for j=1:length(path)-1
            count = count + 1;
            if count == 1
                swc_tmp = [path(j+1) ep(seed_idx)];
            else
                tmp = [path(j+1) path(j)];
                [C,~,~] = intersect(tmp, swc_tmp,'rows');
                if isempty(C)
                    swc_tmp = [swc_tmp; tmp];
                end
            end
        end
    end


    %!!!Don't use unique, it will sort you data, and will lose spatial information.
    % swc_tmp = unique(swc_tmp, 'rows');

    swc_tmp = [list_neuron(swc_tmp(:, 1)), list_neuron(swc_tmp(:, 2))];
    % add root idx
    swc_tmp = [[list_neuron(ep(seed_idx)) -1]; swc_tmp];

    idx_tmp = zeros(size(swc_tmp,1),1);

    for i=1:size(swc_tmp,1)
        if swc_tmp(i,2)==-1
            idx_tmp(i) = -1;
        else
            parent_idx = find(swc_tmp(:,1)==swc_tmp(i,2), 1);
            if isempty(parent_idx)
                swc_tmp(i,2)
            end
            idx_tmp(i) = parent_idx;
        end
    end
    swc_tmp = [swc_tmp(:,1) swc_tmp(:,2) idx_tmp];
    [x,y,z] = ind2sub([w+2,l+2,h+2], swc_tmp(:,1));
    swc = [(1:size(x,1))', zeros(size(x,1),1), x-1, y-1, z-1, zeros(size(x,1),1), swc_tmp(:,3)];
    
    [x,y,z]=ind2sub([w+2,l+2,h+2], list_neuron(ep(seed_idx)));
    seed_xyz=[x-1,y-1,z-1];
else
    swc=[];
    seed_xyz=[];
end
end

