function run_on_simData(exp_name, config)
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
    exp_name='STRCF_sim_data/';
end

sim_dir='/path/to/xBTracer/simulation/simulated_data/';

tif_dir='/path/to/xBTracer/simulation/simulated_data/*.tif';
tifs=dir(tif_dir);
tifs=natsortfiles(tifs);

for tif_id=1:numel(tifs)
    tif_name=split(tifs(tif_id).name, 'simVol_');
    tif_name=tif_name(2);
    tif_name=split(tif_name, '.tif');
    tif_name=tif_name(1);
    tif_name=tif_name{1};
    
    % read bbs
    bb_dir = [sim_dir 'bbs/' tif_name '.txt'];
    bb_cells = read_bbs(bb_dir);

    % save_dir
    result_dir = '/path/to//xBTracer/results/';
    result_dir = [result_dir exp_name tif_name '/'];

    [imgs,h,w] = imread_big(fullfile(tifs(tif_id).folder,tifs(tif_id).name));

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
            if bb_cells{idx,3}(3)<=2 || bb_cells{idx,3}(4)<=2 || prod(bb_cells{idx,3}(3:4)) >= 25*25
                continue;
            end

    %         init_params = read_default_csr_parameters();
    %         writestruct(init_params, config.param_file)

            swc=single_tracing(bb_cells{idx,3}, bb_cells{idx,2}, imgs, config);
            if ~isempty(swc)
                swc(1,end)=-1;
            end
            save_to_swc_file(swc, config.swc_file);
        end

        % eval
        score_struct=readstruct(xml_file);
        
        target_swc_dir=strcat(sim_dir, 'swc_gt/', tif_name, '/', bb_cells{idx,1}, '.swc');
        test_swc_dir=strcat(result_dir, 'swc/', bb_cells{idx,1}, '.swc');
        
        if config.visualize_tracker
            plot_swcs({target_swc_dir, test_swc_dir});
        end

        score_struct=save_eval_xml(score_struct, bb_cells{idx,1}, test_swc_dir, target_swc_dir);

        % write back
        final_strcut.metrics=score_struct;
        struct2xml(final_strcut, xml_file);
    end
end