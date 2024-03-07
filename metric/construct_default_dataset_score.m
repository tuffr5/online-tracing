function construct_default_dataset_score(bb_dir, xml_file)
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

cnt=1;
for i=1:length(bblist)
    tmp=split(bblist{i}, ' ');
    if ~strcmp(tmp(3), '-')
        id=strcat('n_', num2str(cnt));
        params.metrics.(id).name=tmp{1};
        params.metrics.(id).vsa=0;
        params.metrics.(id).MFD=0;
        params.metrics.(id).BRR=0;
        params.metrics.(id).length=0;
        params.metrics.(id).num_bif=0;
        cnt=cnt+1;
    end
end

struct2xml(params, xml_file);
end