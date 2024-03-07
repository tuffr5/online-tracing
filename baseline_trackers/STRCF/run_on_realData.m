function run_on_realData(exp_name, config)
PLANE=['xy' 'zy' 'xz']; % z-1, x-2, y-3

% set this to tracker directory
tracker_path='./';
% add paths
addpath(tracker_path);
addpath(fullfile(tracker_path, 'utils'));

setup_paths();

config.enable_bifurcation=false;
config.enable_cross_plane=false;
config.enable_alpha_mask=false;
config.visualize_tracker=false;
config.debug=false;

if ~config.enable_bifurcation
    config.visualize_bifurcation=false;
end

if nargin < 1
    exp_name='STRCF_real_data/';
end

% read bbs
bb_dir = '/path/to/xBTracer/dataset/gt_bb.txt';
bb_cells = read_bbs(bb_dir);

% save_dir
result_dir = '/path/to/xBTracer/results/';
result_dir = [result_dir exp_name];

[imgs,h,w] = imread_big('/path/to/nTracer_sample_denoised.tif');

%Height*Width*(CZ) --> Height*Width*C*Z !!!important
imgs = reshape(imgs,h,w,3,[]);
% YXCZ -> ZYXC
imgs=permute(imgs,[4 1 2 3]);

[d,h,w,c]=size(imgs);

for idx=1:length(bb_cells)
    config.txt_file=[result_dir 'raw/' bb_cells{idx,1} '.txt'];
    config.swc_file=[result_dir 'swc/' bb_cells{idx,1} '.swc'];
    
    xml_file=[result_dir 'score.xml'];
    
    if idx==1            
        construct_default_dataset_score(bb_dir, xml_file);
    end

    if ~isfile(config.swc_file)
        swc=single_tracing(bb_cells{idx,3}, bb_cells{idx,2}, imgs, config);
        if ~isempty(swc)
            swc(1,end)=-1;
        end
        save_to_swc_file(swc, config.swc_file);
    end
    
    % construct_default_dataset_score(bb_dir, xml_file);
    score_struct=readstruct(xml_file);

    target_swc_dir=strcat('/path/to/xBTracer/dataset/gt/', bb_cells{idx,1}, '.swc');
    test_swc_dir=strcat(result_dir, 'swc/', bb_cells{idx,1}, '.swc');

    score_struct=save_eval_xml(score_struct, bb_cells{idx,1}, test_swc_dir, target_swc_dir);

    % write back
    final_strcut.metrics=score_struct;
    struct2xml(final_strcut, xml_file);
end
end