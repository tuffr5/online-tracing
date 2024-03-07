bb_dir='simulated_data/bak/*.txt';
result_dir='../results/xBTracer_sim_data_XYZ/';
new_result_dir='../results/xBTracer_sim_data/';

newbb_dir='simulated_data/bbs/';

bbFiles=dir(bb_dir);


for bbid=1:numel(bbFiles)
	sub_dir=split(bbFiles(bbid).name, '.txt');
	sub_dir=sub_dir(1);
	sub_dir=sub_dir{1};
    
    num_cells=split(sub_dir, '_');
    num_cells=num_cells(1);
    num_cells=split(num_cells{1}, 'Cell');
    num_cells=num_cells(2);
    num_cells=str2num(num_cells{1});

	% new bbs
	fid=fopen(fullfile(newbb_dir, bbFiles(bbid).name), 'w');

	result_xml=fullfile(result_dir, sub_dir, 'score.xml');

	% load bbs
	bb_cells = read_bbs(fullfile(bbFiles(bbid).folder, bbFiles(bbid).name));

	% LOAD scores
	score_struct=readstruct(result_xml);

	for j=1:num_cells
		% get the scores
		scores=[0 0 0];

		fn = fieldnames(score_struct);
		for i=1:numel(fn)
		    if( strcmp(score_struct.(fn{i}).name, [num2str(j) '-xy']) )
		        scores(1)=score_struct.(fn{i}).vsa;
		    end

		    if( strcmp(score_struct.(fn{i}).name, [num2str(j) '-zy']) )
		        scores(2)=score_struct.(fn{i}).vsa;
		    end

		    if( strcmp(score_struct.(fn{i}).name, [num2str(j) '-xz']) )
		        scores(3)=score_struct.(fn{i}).vsa;
		    end
		end
		% get the biggest one

        if max(scores(:)) == 0
        	continue;
        end
        
		ids=find(scores==max(scores(:)));
		id=ids(1);

		switch id
            case 1
                real_id_name=[num2str(j) '-xy'];
            case 2
                real_id_name=[num2str(j) '-zy'];
            case 3
                real_id_name=[num2str(j) '-xz'];
            otherwise
                error('No sich id.');
        end
        
        for ii=1:numel(bb_cells)
            if strcmp(bb_cells{ii, 1}, real_id_name)
                real_id=ii;
                break;
            end
        end

		% write the corresponding bb to newbb_file
		thisbb=sprintf('%s %s %d %d %d %d %d %d\n', num2str(j), bb_cells{real_id,2}, bb_cells{real_id,3});
		fprintf(fid, thisbb);

		% move and rename swc file
		swc_file=fullfile(result_dir, sub_dir, 'swc', [bb_cells{real_id,1}, '.swc']);
		new_swc_file=fullfile(new_result_dir, sub_dir, 'swc', [num2str(j), '.swc']);

		copyfile(swc_file, new_swc_file);

		txt_file=fullfile(result_dir, sub_dir, 'raw', [bb_cells{real_id,1}, '.txt']);
		new_txt_file=fullfile(new_result_dir, sub_dir, 'raw', [num2str(j), '.txt']);

		copyfile(txt_file, new_txt_file);
    end
    
    fclose(fid);
end