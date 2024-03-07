function [tracker, region, curr_score]=track_STRCF_tracker(tracker, im, frame)
% pack up tracker
pos=tracker.pos;
params=tracker.params;
currentScaleFactor=tracker.currentScaleFactor;
scaleFactors=tracker.scaleFactors;
features=tracker.features;
global_fparams=tracker.global_fparams;
feature_extract_info=tracker.feature_extract_info;
cos_window=tracker.cos_window;

k1=tracker.k1;
cf_f=tracker.cf_f;
f_pre_f=tracker.f_pre_f;
block_inds=tracker.block_inds;
output_sz=tracker.output_sz;
ky=tracker.ky;
kx=tracker.kx;
img_support_sz=tracker.img_support_sz;
min_scale_factor=tracker.min_scale_factor;
max_scale_factor=tracker.max_scale_factor;
yf=tracker.yf;
reg_window=tracker.reg_window;
base_target_sz=tracker.base_target_sz;

old_pos = inf(size(pos));
iter = 1;

% Load learning parameters
admm_max_iterations = params.max_iterations;
init_penalty_factor = params.init_penalty_factor;
max_penalty_factor = params.max_penalty_factor;
penalty_scale_step = params.penalty_scale_step;
temporal_regularization_factor = params.temporal_regularization_factor; 
newton_iterations = params.newton_iterations;

%translation search
while iter <= params.refinement_iterations && any(old_pos ~= pos)
    % Extract features at multiple resolutions
    sample_pos = round(pos);
    sample_scale = currentScaleFactor*scaleFactors;
    xt = extract_features(im, sample_pos, sample_scale, features, global_fparams, feature_extract_info);
                            
    % Do windowing of features
    xtw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xt, cos_window, 'uniformoutput', false);
    
    % Compute the fourier series
    xtf = cellfun(@fft2, xtw, 'uniformoutput', false);
                
    % Compute convolution for each feature block in the Fourier domain
    % and the sum over all blocks.
    scores_fs_feat{k1} = gather(sum(bsxfun(@times, conj(cf_f{k1}), xtf{k1}), 3));
    scores_fs_sum = scores_fs_feat{k1};
    for k = block_inds
        scores_fs_feat{k} = gather(sum(bsxfun(@times, conj(cf_f{k}), xtf{k}), 3));
        scores_fs_feat{k} = resizeDFT2(scores_fs_feat{k}, output_sz);
        scores_fs_sum = scores_fs_sum +  scores_fs_feat{k};
    end
     
    % Also sum over all feature blocks.
    % Gives the fourier coefficients of the convolution response.
    scores_fs = permute(gather(scores_fs_sum), [1 2 4 3]);
    
    responsef_padded = resizeDFT2(scores_fs, output_sz);
    response = ifft2(responsef_padded, 'symmetric');
    
    % test target presentence for occlusion
    curr_score = [];
    resp_quality = max(response(:));
    
    if isempty(tracker.resp_budg)
        % in first localization frame response score needs to be added as
        % response normalization score
        % therefore 1 is added into the budget
        tracker.resp_norm = resp_quality;
        tracker.resp_budg(end+1) = 1;
    else
        % normalize current response score and add it to the budget
        tracker.resp_budg(end+1) = resp_quality / tracker.resp_norm;
        % if budget is reached, remove first (the oldest) element
        if numel(tracker.resp_budg) > tracker.resp_budg_sz
            tracker.resp_budg(1) = [];
        end
    end
    
    if numel(tracker.resp_budg) == tracker.resp_budg_sz
        response_budget_mean = mean(tracker.resp_budg);
        curr_quality_norm = resp_quality / tracker.resp_norm;
        curr_score = (response_budget_mean - curr_quality_norm) / curr_quality_norm;
        if curr_score > tracker.detect_failure
            region=[];  % target lost
            return;
        end
    end
    
    [disp_row, disp_col, sind] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, output_sz);
                
    % Compute the translation vector in pixel-coordinates and round
    % to the closest integer pixel.
    translation_vec = [disp_row, disp_col] .* (img_support_sz./output_sz) * currentScaleFactor * scaleFactors(sind);            
    scale_change_factor = scaleFactors(sind);
    
    % update position
    old_pos = pos;
    pos = sample_pos + translation_vec;
    
    if params.clamp_position
        pos = max([1 1], min([size(im,1) size(im,2)], pos));
    end
                
    % Update the scale
    currentScaleFactor = currentScaleFactor * scale_change_factor;
    
    % Adjust to make sure we are not to large or to small
    if currentScaleFactor < min_scale_factor
        currentScaleFactor = min_scale_factor;
    elseif currentScaleFactor > max_scale_factor
        currentScaleFactor = max_scale_factor;
    end
    
    iter = iter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model update step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract image region for training sample
if isnan(pos(1)) || isnan(pos(2))
    region=[];
    return;
end
sample_pos = round(pos);
xl = extract_features(im, sample_pos, currentScaleFactor, features, global_fparams, feature_extract_info);

% do windowing of features
xlw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xl, cos_window, 'uniformoutput', false);

% compute the fourier series
xlf = cellfun(@fft2, xlw, 'uniformoutput', false);

% train the CF model for each feature
for k = 1: numel(xlf)
    model_xf = xlf{k};

    mu = temporal_regularization_factor(k);
    
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
sz = tracking_result.target_size;
tl = tracking_result.center_pos - (sz - 1)/2;
br = tracking_result.center_pos + (sz - 1)/2;
x1 = tl(2); y1 = tl(1);
x2 = br(2); y2 = br(1);

region=round(double([x1, y1, x2-x1+1, y2-y1+1]));

% set to tracker
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
end