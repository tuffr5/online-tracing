function bb_cells = read_bbs(bb_dir)
bblist = [];
fid = fopen(bb_dir);
if fid==-1
    disp(['Error to open the file : ' filename]);
    return;
else
    i=1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        bblist{i} = deblank(tline);
        i = i+1;
    end
end
fclose(fid);

bb_cells = {};
cnt=1;
for i=1:length(bblist)
    tmp=split(bblist{i}, ' ');
    if ~strcmp(tmp(3), '-') && ~strcmp(tmp(4), '-')
        bb_cells{cnt,1}=tmp{1}; % neuron_id
        bb_cells{cnt,2}=tmp{2}; % direction
        bb_cells{cnt,3}=[str2num(tmp{3}) str2num(tmp{4}) str2num(tmp{5}) ...
            str2num(tmp{6}) str2num(tmp{7}) str2num(tmp{8})]; % bb6
        cnt=cnt+1;
    end
end

end

