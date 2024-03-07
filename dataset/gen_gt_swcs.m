clear all;
clc;

filename='/home/path/to/Downloads/nTracer_sample_img/nTracer sample.tif-data.txt';
out_dir='~/Downloads/experiments/csr-dcf/dataset/gt/';

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
    
    if findstr('Spine', tline(1:5))
        break;
    end 

%     if ~isempty(neurite) && ~startsWith(neurite, 'Neurite 5-3')
%         continue;
%     end
    
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

%         data{i,1} = neuron;
%         if findstr('Soma', point(idx_p(1) + 1:idx_p(2) - 1))
%             data{i,2} = soma;
%         else
%             data{i,2} = neurite;
%         end
        
        % only save neurite
        
        if findstr('Soma', point(idx_p(1) + 1:idx_p(2) - 1))
            continue;
        else
            tmp_cell=split(neurite, ' ');
            data{i,1} = tmp_cell{2};
            data{i,2} = p;
            i = i + 1;
        end
    end
end
fclose(fid);

ids=[];

for i=1:size(data,1)
    if length(split(data{i,1}, '-'))==2
        ids=[ids; i];
    end
end

for i=1:size(ids,1)-1
    p=vertcat(data(ids(i):ids(i+1)-1,2));
    save_gt_tracing(cell2mat(p), [out_dir data{ids(i), 1} '.swc']);
end

p=vertcat(data(ids(end):end,2));
save_gt_tracing(cell2mat(p), [out_dir data{ids(end), 1} '.swc']);