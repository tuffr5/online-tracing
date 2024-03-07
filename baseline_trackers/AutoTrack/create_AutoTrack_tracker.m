function tracker=create_AutoTrack_tracker(im, bb)
params=read_default_params();
frame=1;


% since hog cell_size is 4 minimum
if bb(4) < 4
    bb(4) = 4;
end

if bb(3) < 4
    bb(3) = 4;
end


params.wsize    = [bb(4), bb(3)];
params.init_pos = [bb(2), bb(1)] + floor(params.wsize/2);

%   Setting parameters for local use.
admm_iterations = params.admm_iterations;
search_area_scale   = params.search_area_scale;
max_image_sample_size=params.max_image_sample_size;
min_image_sample_size=params.min_image_sample_size;
output_sigma_factor = params.output_sigma_factor;
% Scale parameters
num_scales=params.num_scales;
scale_sigma_factor=params.scale_sigma_factor;
scale_step=params.scale_step;
scale_lambda=params.scale_lambda;   
scale_model_factor=params.scale_model_factor;
scale_model_max_area =params.scale_model_max_area;
lambda=params.admm_lambda;
features    = params.t_features;
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);
visualization  = params.visualization;

epsilon=params.epsilon;
delta=params.delta;
zeta=params.zeta;
newton_iterations = params.newton_iterations;
featureRatio = params.t_global.cell_size;
search_area = prod(target_sz * search_area_scale);
global_feat_params = params.t_global;
nu=params.nu;
if search_area > max_image_sample_size
    currentScaleFactor = sqrt(search_area / max_image_sample_size);
elseif search_area <min_image_sample_size
    currentScaleFactor = sqrt(search_area / min_image_sample_size);
else
    currentScaleFactor = 1.0;
end
% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;
reg_sz= floor(base_target_sz/featureRatio);
% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor(base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
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
[rs, cs]     = ndgrid(rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y); %   FFT of y.\

interp_sz = use_sz;

% construct cosine window
cos_window = single(hann(use_sz(1)+2)*hann(use_sz(2)+2)');
cos_window = cos_window(2:end-1,2:end-1);

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

%% SCALE ADAPTATION INITIALIZATION
% Use the translation filter to estimate the scale
scale_sigma = sqrt(num_scales) * scale_sigma_factor;
ss = (1:num_scales) - ceil(num_scales/2);
ys = exp(-0.5 * (ss.^2) / scale_sigma^2);
ysf = single(fft(ys));
if mod(num_scales,2) == 0
    scale_window = single(hann(num_scales+1));
    scale_window = scale_window(2:end);
else
    scale_window = single(hann(num_scales));
end
ss = 1:num_scales;
scaleFactors = scale_step.^(ceil(num_scales/2) - ss);
if scale_model_factor^2 * prod(target_sz) > scale_model_max_area
    scale_model_factor = sqrt(scale_model_max_area/prod(target_sz));
end
if prod(target_sz) >scale_model_max_area
    params.scale_model_factor = sqrt(scale_model_max_area/prod(target_sz));
end
scale_model_sz = floor(target_sz * scale_model_factor);

% set maximum and minimum scales
min_scale_factor = scale_step ^ ceil(log(max(5 ./sz)) / log(scale_step));
max_scale_factor =scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
    
% Pre-computes the grid that is used for score optimization
ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';

if size(im,3) > 1 && colorImage == false
    im = im(:,:,1);
end

if frame==1   
% extract training sample image region
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    pixels = uint8(gather(pixels));
    x=get_features(pixels,features,global_feat_params);
    xf=fft2(bsxfun(@times,x,cos_window));
else
% use detection features
     shift_samp_pos = 2*pi * translation_vec ./(currentScaleFactor* sz);
     xf = shift_sample(xtf, shift_samp_pos, kx', ky');
end
        
if  frame == 1
    [range_h,range_w,w]=init_regwindow(use_sz,reg_sz,params);
    g_pre= zeros(size(xf));
    mu = 0;
else
    mu=zeta;
end

occ=false;

if ~occ
    g_f = single(zeros(size(xf)));
    h_f = g_f;
    l_f = h_f;
    gamma = 1;
    betha = 10;
    gamma_max = 10000;
    
    
    % ADMM solution    
    T = prod(use_sz);
    S_xx = sum(conj(xf) .* xf, 3);
    Sg_pre= sum(conj(xf) .* g_pre, 3);
    Sgx_pre= bsxfun(@times, xf, Sg_pre);
    iter = 1;
    while (iter <= admm_iterations)
        % subproblem g
        B = S_xx + T * (gamma + mu);
        Shx_f = sum(conj(xf) .* h_f, 3);
        Slx_f = sum(conj(xf) .* l_f, 3);
        g_f = ((1/(T*(gamma + mu)) * bsxfun(@times,  yf, xf)) - ((1/(gamma + mu)) * l_f) +(gamma/(gamma + mu)) * h_f) + (mu/(gamma + mu)) * g_pre - ...
            bsxfun(@rdivide,(1/(T*(gamma + mu)) * bsxfun(@times, xf, (S_xx .*  yf)) + (mu/(gamma + mu)) * Sgx_pre- ...
            (1/(gamma + mu))* (bsxfun(@times, xf, Slx_f)) +(gamma/(gamma + mu))* (bsxfun(@times, xf, Shx_f))), B);
        %   subproblem h
        lhd= T ./  (lambda*w .^2 + gamma*T); 
        X=ifft2(gamma*(g_f + l_f));
        h=bsxfun(@times,lhd,X);
        h_f = fft2(h);
        %   subproblem mu
        if frame>2&&iter<admm_iterations
            for i=1:size(g_f,3)
                z=power(norm(g_f(:,:,i)-g_pre(:,:,i),2),2)/(2*epsilon);
                mu=ref_mu-z;
            end
        end
        %   update h
        l_f = l_f + (gamma * (g_f - h_f));
        %   update gamma
        gamma = min(betha* gamma, gamma_max);
        iter = iter+1;
     end
            
end
% save the trained filters
g_pre= g_f;

%% Upadate Scale
if frame==1
    xs = crop_scale_sample(im, pos, base_target_sz, currentScaleFactor * scaleFactors, scale_window, scale_model_sz);
else
    xs= shift_sample_scale(im, pos, base_target_sz,xs,recovered_scale,currentScaleFactor*scaleFactors,scale_window,scale_model_sz);
end

xsf = fft(xs,[],2);
new_sf_num = bsxfun(@times, ysf, conj(xsf));
new_sf_den = sum(xsf .* conj(xsf), 1);

if frame == 1
    sf_den = new_sf_den;
    sf_num = new_sf_num;
else
    sf_den = (1 - params.learning_rate_scale) * sf_den + params.learning_rate_scale * new_sf_den;
    sf_num = (1 - params.learning_rate_scale) * sf_num + params.learning_rate_scale * new_sf_num;
end

tracker.params=params;
tracker.c=[params.init_pos(2) params.init_pos(1)];
tracker.pos=params.init_pos;


tracker.g_f=g_f;
tracker.sf_num=sf_num;
tracker.sf_den=sf_den;
tracker.g_pre=g_pre;
tracker.w=w;
tracker.range_h=range_h;
tracker.range_w=range_w;

% target presentence
tracker.resp_budg_sz = params.resp_budg_sz;
tracker.resp_budg = [];
tracker.resp_norm = 0;
tracker.detect_failure = params.detect_failure;
tracker.cutoff_intensity=params.cutoff_intensity;
end