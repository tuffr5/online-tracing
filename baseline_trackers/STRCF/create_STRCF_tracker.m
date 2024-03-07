function tracker=create_STRCF_tracker(im, bb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params=read_default_params();

cx=bb(1)+0.5*bb(3);
cy=bb(2)+0.5*bb(4);

% Init position
pos = [cy cx];
target_sz = [bb(4) bb(3)];
params.init_sz = target_sz;

% Feature settings
features = params.t_features;

% Set default parameters
params = init_default_params(params);

% Global feature parameters
if isfield(params, 't_global')
    global_fparams = params.t_global;
else
    global_fparams = [];
end

global_fparams.use_gpu = params.use_gpu;
global_fparams.gpu_id = params.gpu_id;

% Define data types
if params.use_gpu
    params.data_type = zeros(1, 'single', 'gpuArray');
else
    params.data_type = zeros(1, 'single');
end
params.data_type_complex = complex(params.data_type);

global_fparams.data_type = params.data_type;

% Load learning parameters
admm_max_iterations = params.max_iterations;
init_penalty_factor = params.init_penalty_factor;
max_penalty_factor = params.max_penalty_factor;
penalty_scale_step = params.penalty_scale_step;
temporal_regularization_factor = params.temporal_regularization_factor; 

init_target_sz = target_sz;

% Check if color image
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        is_color_image = false;
    else
        is_color_image = true;
    end
else
    is_color_image = false;
end

if size(im,3) > 1 && is_color_image == false
    im = im(:,:,1);
end

% Check if mexResize is available and show warning otherwise.
params.use_mexResize = true;
global_fparams.use_mexResize = true;
try
    [~] = mexResize(ones(5,5,3,'uint8'), [3 3], 'auto');
catch err
    params.use_mexResize = false;
    global_fparams.use_mexResize = false;
end

% Calculate search area and initial scale factor
search_area = prod(init_target_sz * params.search_area_scale);
if search_area > params.max_image_sample_size
    currentScaleFactor = sqrt(search_area / params.max_image_sample_size);
elseif search_area < params.min_image_sample_size
    currentScaleFactor = sqrt(search_area / params.min_image_sample_size);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        img_sample_sz = floor(base_target_sz * params.search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        img_sample_sz = repmat(sqrt(prod(base_target_sz * params.search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        img_sample_sz = base_target_sz + sqrt(prod(base_target_sz * params.search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    case 'custom'
        img_sample_sz = [base_target_sz(1)*2 base_target_sz(2)*2];
end

[features, global_fparams, feature_info] = init_features(features, global_fparams, is_color_image, img_sample_sz, 'exact');

% Set feature info
img_support_sz = feature_info.img_support_sz;
feature_sz = unique(feature_info.data_sz, 'rows', 'stable');
feature_cell_sz = unique(feature_info.min_cell_size, 'rows', 'stable');
num_feature_blocks = size(feature_sz, 1);

% Get feature specific parameters
feature_extract_info = get_feature_extract_info(features);

% Size of the extracted feature maps
feature_sz_cell = mat2cell(feature_sz, ones(1,num_feature_blocks), 2);
filter_sz = feature_sz;
filter_sz_cell = permute(mat2cell(filter_sz, ones(1,num_feature_blocks), 2), [2 3 1]);

% The size of the label function DFT. Equal to the maximum filter size
[output_sz, k1] = max(filter_sz, [], 1);
k1 = k1(1);

% Get the remaining block indices
block_inds = 1:num_feature_blocks;
block_inds(k1) = [];

% Construct the Gaussian label function
yf = cell(numel(num_feature_blocks), 1);
for i = 1:num_feature_blocks
    sz = filter_sz_cell{i};
    output_sigma = sqrt(prod(floor(base_target_sz/feature_cell_sz(i)))) * params.output_sigma_factor;
    rg           = circshift(-floor((sz(1)-1)/2):ceil((sz(1)-1)/2), [0 -floor((sz(1)-1)/2)]);
    cg           = circshift(-floor((sz(2)-1)/2):ceil((sz(2)-1)/2), [0 -floor((sz(2)-1)/2)]);
    [rs, cs]     = ndgrid(rg,cg);
    y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
    yf{i}           = fft2(y); 
end

% Compute the cosine windows
cos_window = cellfun(@(sz) hann(sz(1))*hann(sz(2))', feature_sz_cell, 'uniformoutput', false);

% Define spatial regularization windows
reg_window = cell(num_feature_blocks, 1);
for i = 1:num_feature_blocks
    reg_scale = floor(base_target_sz/params.feature_downsample_ratio(i));
    use_sz = filter_sz_cell{i};    
    reg_window{i} = ones(use_sz) * params.reg_window_max;
    range = zeros(numel(reg_scale), 2);
    
    % determine the target center and range in the regularization windows
    for j = 1:numel(reg_scale)
        range(j,:) = [0, reg_scale(j) - 1] - floor(reg_scale(j) / 2);
    end
    center = floor((use_sz + 1)/ 2) + mod(use_sz + 1,2);
    range_h = (center(1)+ range(1,1)) : (center(1) + range(1,2));
    range_w = (center(2)+ range(2,1)) : (center(2) + range(2,2));
    
    % some bugs for Brainbow images since target is too small
    range_h(range_h<1)=[];
    range_h(range_h>size(reg_window{i},1))=[];
    
    range_w(range_w<1)=[];
    range_w(range_w>size(reg_window{i},2))=[];
    
    reg_window{i}(range_h, range_w) = params.reg_window_min;
end

% Pre-computes the grid that is used for socre optimization
ky = circshift(-floor((filter_sz_cell{1}(1) - 1)/2) : ceil((filter_sz_cell{1}(1) - 1)/2), [1, -floor((filter_sz_cell{1}(1) - 1)/2)]);
kx = circshift(-floor((filter_sz_cell{1}(2) - 1)/2) : ceil((filter_sz_cell{1}(2) - 1)/2), [1, -floor((filter_sz_cell{1}(2) - 1)/2)])';
newton_iterations = params.newton_iterations;

% Use the translation filter to estimate the scale
nScales = params.number_of_scales;
scale_step = params.scale_step;
scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
scaleFactors = scale_step .^ scale_exp;

if nScales > 0
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ img_support_sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

seq.time = 0;

% Define the learning variables
f_pre_f = cell(num_feature_blocks, 1);
cf_f = cell(num_feature_blocks, 1);

% Allocate
scores_fs_feat = cell(1,1,num_feature_blocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model update step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract image region for training sample
sample_pos = round(pos);
xl = extract_features(im, sample_pos, currentScaleFactor, features, global_fparams, feature_extract_info);

% do windowing of features
xlw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xl, cos_window, 'uniformoutput', false);

% compute the fourier series
xlf = cellfun(@fft2, xlw, 'uniformoutput', false);

% train the CF model for each feature
for k = 1: numel(xlf)
    model_xf = xlf{k};

    f_pre_f{k} = zeros(size(model_xf));
    mu = 0;


    % intialize the variables
    f_f = single(zeros(size(model_xf)));
    g_f = f_f;
    h_f = f_f;
    gamma  = init_penalty_factor(k);
    gamma_max = max_penalty_factor(k);
    gamma_scale_step = penalty_scale_step(k);

    % use the GPU mode
    if params.use_gpu
        model_xf = gpuArray(model_xf);
        f_f = gpuArray(f_f);
        f_pre_f{k} = gpuArray(f_pre_f{k});
        g_f = gpuArray(g_f);
        h_f = gpuArray(h_f);
        reg_window{k} = gpuArray(reg_window{k});
        yf{k} = gpuArray(yf{k});
    end

    % pre-compute the variables
    T = prod(output_sz);
    S_xx = sum(conj(model_xf) .* model_xf, 3);
    Sf_pre_f = sum(conj(model_xf) .* f_pre_f{k}, 3);
    Sfx_pre_f = bsxfun(@times, model_xf, Sf_pre_f);

    % solve via ADMM algorithm
    iter = 1;
    while (iter <= admm_max_iterations)

        % subproblem f
        B = S_xx + T * (gamma + mu);
        Sgx_f = sum(conj(model_xf) .* g_f, 3);
        Shx_f = sum(conj(model_xf) .* h_f, 3);

        f_f = ((1/(T*(gamma + mu)) * bsxfun(@times,  yf{k}, model_xf)) - ((1/(gamma + mu)) * h_f) +(gamma/(gamma + mu)) * g_f) + (mu/(gamma + mu)) * f_pre_f{k} - ...
            bsxfun(@rdivide,(1/(T*(gamma + mu)) * bsxfun(@times, model_xf, (S_xx .*  yf{k})) + (mu/(gamma + mu)) * Sfx_pre_f - ...
            (1/(gamma + mu))* (bsxfun(@times, model_xf, Shx_f)) +(gamma/(gamma + mu))* (bsxfun(@times, model_xf, Sgx_f))), B);

        %   subproblem g
        g_f = fft2(argmin_g(reg_window{k}, gamma, real(ifft2(gamma * f_f+ h_f)), g_f));

        %   update h
        h_f = h_f + (gamma * (f_f - g_f));

        %   update gamma
        gamma = min(gamma_scale_step * gamma, gamma_max);

        iter = iter+1;
    end

    % save the trained filters
    f_pre_f{k} = f_f;
    cf_f{k} = f_f;
end  

% Update the target size (only used for computing output box)
target_sz = base_target_sz * currentScaleFactor;

%save position and calculate FPS
tracking_result.center_pos = double(pos);
tracking_result.target_size = double(target_sz);

% pack up to tracker
tracker.pos=pos;
tracker.params=params;
tracker.currentScaleFactor=currentScaleFactor;
tracker.scaleFactors=scaleFactors;
tracker.features=features;
tracker.global_fparams=global_fparams;
tracker.feature_extract_info=feature_extract_info;
tracker.cos_window=cos_window;

tracker.k1=k1;
tracker.cf_f=cf_f;
tracker.f_pre_f=f_pre_f;
tracker.c=[tracking_result.center_pos(2) tracking_result.center_pos(1)];
tracker.block_inds=block_inds;
tracker.output_sz=output_sz;
tracker.ky=ky;
tracker.kx=kx;

tracker.img_support_sz=img_support_sz;
tracker.min_scale_factor=min_scale_factor;
tracker.max_scale_factor=max_scale_factor;
tracker.yf=yf;
tracker.reg_window=reg_window;
tracker.base_target_sz=base_target_sz;

% target presentence
tracker.resp_budg_sz = params.resp_budg_sz;
tracker.resp_budg = [];
tracker.resp_norm = 0;
tracker.detect_failure = params.detect_failure;
tracker.cutoff_intensity=params.cutoff_intensity;
end