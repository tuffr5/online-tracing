function parameters = read_default_csr_parameters(p)

    parameters.save_output_stats = false;

    % filter parameters
    parameters.padding = 2;
    parameters.learning_rate = 0.01;
%     parameters.feature_type = {'hog', 'cn', 'gray'};
    parameters.feature_type = {'cn', 'gray'}; % only use h in hsv
    parameters.y_sigma = 1;
    parameters.channels_weight_lr = parameters.learning_rate;
    parameters.use_channel_weights = true;
    
    % target presentence
    parameters.resp_budg_sz = 20;
    parameters.detect_failure = 3; % to check object existance
    
    % mask history
    parameters.mask_budg_sz = 50;
    
    % cross section reliablity
    parameters.cross_plane_degree = 0.9; % lower bound to change plane
    parameters.maskbb_size_sigma = 0.4;
    parameters.maskbb_hwratio_sigma = 0.4;
    
    % for branch backtracking path cost
    parameters.branch_path_degree = 0.9; % lower bound to add branch
    parameters.branch_color_sigma = 0.2;
    parameters.branch_bb_padding = 2; % compensation for unsaturated mask
    parameters.branch_counter = 3; % how many times appeared to count as a branch

    % segmentation parameters
    parameters.hist_lr = 0.01;
    parameters.retain_rate = 0.2;
    parameters.nbins = 32;  % N bins for segmentation
    parameters.seg_colorspace = 'rgb';     % 'rgb' or 'hsv'
    parameters.use_segmentation = true;  % false to disable use of segmentation
    parameters.mask_diletation_type = 'disk';  % for function strel (square, disk, ...)
    parameters.mask_diletation_sz = 1;
    parameters.cutoff_intensity = 15; % exclude low intensity pixels
    parameters.enable_alpha_mask=true;

    % scale adaptation parameters (from DSST)
    parameters.currentScaleFactor = 1.0;
    parameters.n_scales = 30;
    parameters.scale_model_factor = 1.0;
    parameters.scale_sigma_factor = 1/4;
    parameters.scale_step = 1.03;
    parameters.scale_model_max_area = 32*16;
    parameters.scale_lr = 0.2;
    
    % diameter adaptation parameters

    % overwrite parameters that come frome input argument
    if nargin > 0
    	fields = fieldnames(p);

    	for i=1:numel(fields)
            if ~isfield(parameters, fields{i})
                warning('Setting parameter value for: %s. It is not set by default.', fields{i});
            end
    		parameters = setfield(parameters, fields{i}, p.(fields{i}));
    	end
    end

end  % endfunction
