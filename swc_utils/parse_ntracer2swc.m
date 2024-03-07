function parse_ntracer2swc(varargin)
%--------------------------------------------------------------------------
% Date: July 2, 2020
% input : tracing file (.txt) generated from nTracer - ImageJ plugin
% output: cell array of input
% % first column : Neuron no.
% % second column: Soma/Neurite
% % third column: xyz coordinates
% organized_data is for 3d matting purpose
%--------------------------------------------------------------------------
if isempty(varargin)
    error("No inpur params");
end

filename=varargin{1};
save_to_swc=varargin{2};
n_id=varargin{3};

fid = fopen(filename);

for k = 1:6
    fgetl(fid);
end

global neuron
global soma
global neurite
% global spine

i = 1;
while not(feof(fid))
    tline = fgetl(fid);

    if tline == -1
        break;
    end

    if findstr('Neuron', tline(1:6))
        num_id = regexp(tline, '\d');
        neuron = str2num(tline(num_id));
    end

    if findstr('Soma', tline(1:4))
        soma = tline;
    end

    if findstr('Neurite', tline(1:7))
        neurite = tline;
    end
    
    if neuron~=n_id
        continue;
    end

    idx = findstr('POINT:', tline);
    if ~isempty(idx)
        p = zeros(length(idx), 3);

        for j = 1:length(idx)
            if j == length(idx)
                point = tline(idx(j):end);
            else
                point = tline(idx(j):idx(j + 1)-1);
            end

            idx_p = findstr(' ', point);
            point_x = str2num(point(idx_p(2) + 1:idx_p(3) - 1));
            point_y = str2num(point(idx_p(3) + 1:idx_p(4) - 1));
            point_z = str2num(point(idx_p(4) + 1:idx_p(5) - 1));

            p(j,:) = [point_x, point_y, point_z];
        end

        data{i,1} = neuron;
        if findstr('Soma', point(idx_p(1) + 1:idx_p(2) - 1))
            data{i,2} = soma;
        else
            data{i,2} = neurite;
        end

        data{i,3} = p;
        i = i + 1;
    end
end
fclose(fid);

tracePlot(data);

ids_pool=[];% no repetetions
points_pool=[]; % no repetetions
pids_pool=[];

id=1;
pid=-1;
for i=1:size(data,1)
    points=data{i,3};
    for j=1:size(points,1)
        % find if in ids_pool
        if ~isempty(points_pool)
            C = intersect(points_pool,points(j,:),'rows');
            if ~isempty(C)
                pid=find(ismember(points_pool, C,'rows'));
                continue;
            end
        end
        
        points_pool=[points_pool; points(j,:)];
        ids_pool=[ids_pool; id];
        pids_pool=[pids_pool; pid];
        id=id+1;
        pid=id-1;
    end
end

swc=[];
for i=1:size(points_pool,1)
    % set to swc
    swc = [swc; [ids_pool(i), 0, points_pool(i,1), points_pool(i,2), points_pool(i,3), 0, pids_pool(i)]];
end

save_to_swc_file(swc, save_to_swc);
end