function tracker=create_ASRCF_tracker(im, bb)
params=read_default_params();
frame=1;

params.wsize    = [bb(4), bb(3)];
params.init_pos = [bb(2), bb(1)] + floor(params.wsize/2);

search_area_scale   = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
learning_rate       = params.learning_rate;
filter_max_area     = params.filter_max_area;
nScales             = params.number_of_scales;
scale_step          = params.scale_step;
interpolate_response = params.interpolate_response;
alphaw=params.alphaw;
update_interval=2;
features    = params.t_features;
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);

visualization  = params.visualization;
init_target_sz = target_sz;
ifcompress=params.ifcompress;
pe=params.pe;

featureRatio = params.t_global.cell_size;
search_area_pos = prod(init_target_sz / featureRatio * search_area_scale);

% when the number of cells are small, choose a smaller cell size
if isfield(params.t_global, 'cell_selection_thresh')
    if search_area_pos < params.t_global.cell_selection_thresh * filter_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));
        
        featureRatio = params.t_global.cell_size;
        search_area_pos = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area_pos > filter_max_area
    currentScaleFactor = sqrt(search_area_pos / filter_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

% construct the label function- correlation output, 2D gaussian function,
% with a peak located upon the target

output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg           = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg           = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs]     = ndgrid( rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y); %   FFT of y.


if interpolate_response == 1
    interp_sz = use_sz * featureRatio;
else
    interp_sz = use_sz;
end

% construct cosine window
feature_sz_cell={use_sz,use_sz,use_sz};
cos_window = cellfun(@(sz) single(hann(sz(1)+2)*hann(sz(2)+2)'), feature_sz_cell, 'uniformoutput', false);
cos_window = cellfun(@(cos_window) cos_window(2:end-1,2:end-1), cos_window, 'uniformoutput', false);

if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

% compute feature dimensionality
feature_dim = 0;
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end
end

if size(im,3) > 1 && colorImage == false
    im = im(:,:,1);
end

if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    scaleFactors = scale_step .^ scale_exp;
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end

% allocate memory for multi-scale tracking
multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');
small_filter_sz = floor(base_target_sz/featureRatio);


if frame==1   
% extract training sample image region
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    pixels = uint8(gather(pixels));
    x=extract_features(pixels,use_sz,features,global_feat_params,frame,ifcompress,pe);
    xf=cellfun(@(feat_map, cos_window) fft2(bsxfun(@times,feat_map,cos_window)), x(1:3), cos_window, 'uniformoutput', false);
    xf=cat(3,xf{1},xf{2},xf{3});
else
% use detection features
    shift_samp_pos = 2*pi * translation_vec ./ (scaleFactors(sind)*currentScaleFactor * sz);
    xf = shift_sample(xtf, shift_samp_pos, kx', ky');
end
xhcf=xf(:,:,1:31);

if (frame == 1)
    model_xf =xf;
    model_xhcf=xhcf;
    model_w=gpuArray(construct_regwindow(use_sz,small_filter_sz));
elseif frame==1||mod(frame,update_interval)==0
    model_xf = ((1 - learning_rate) * model_xf) + (learning_rate * xf);
    % model_xf = xf;
    model_xhcf = ((1 - learning_rate) * model_xhcf) + (learning_rate * xhcf);
end
    
% ADMM solution    
if (frame==1||mod(frame,update_interval)==0) 
    w = gpuArray(params.w_init*single(ones(use_sz)));
% ADMM solution for localization
[g_f,h_f]=ADMM_solve_h(params,use_sz,model_xf,yf,small_filter_sz,w,model_w,frame);
for iteration = 1:params.al_iteration-1
    [w]=ADMM_solve_w(params,use_sz,model_w,h_f);
    [g_f,h_f]=ADMM_solve_h(params,use_sz,model_xf,yf,small_filter_sz,w,model_w,frame);
end
model_w=alphaw*w+(1-alphaw)*model_w;
% ADMM solution for scale estimation
[g_hcf]=ADMM_base(params,use_sz,model_xhcf,xhcf,yf,small_filter_sz,frame);

end

tracker.params=params;
tracker.c=[params.init_pos(2) params.init_pos(1)];
tracker.pos=params.init_pos;
tracker.g_f=g_f;
tracker.g_hcf=g_hcf;
tracker.model_xf=model_xf;
tracker.model_xhcf=model_xhcf;
tracker.yf=yf;
tracker.model_w=model_w;

% target presentence
tracker.resp_budg_sz = params.resp_budg_sz;
tracker.resp_budg = [];
tracker.resp_norm = 0;
tracker.detect_failure = params.detect_failure;
tracker.cutoff_intensity=params.cutoff_intensity;

end