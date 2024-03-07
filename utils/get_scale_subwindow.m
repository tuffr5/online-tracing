function out = get_scale_subwindow(im, pos, base_target_sz, scaleFactors, scale_window, scale_model_sz, scale_feature, w2c)

nScales = length(scaleFactors);

for s = 1:nScales
    patch_sz = floor(base_target_sz * scaleFactors(s));
    
    xs = floor(pos(2)) + (1:patch_sz(2)) - floor(patch_sz(2)/2);
    ys = floor(pos(1)) + (1:patch_sz(1)) - floor(patch_sz(1)/2);
    
    %check for out-of-bounds coordinates, and set them to the values at
    %the borders
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
    
    %extract image
    im_patch = im(ys, xs, :);
    
    % resize image to model size
    im_patch_resized = imresize(im_patch, scale_model_sz, 'bilinear', 'AntiAliasing',false);

    % extract scale features
%     temp_hog = fhog(single(im_patch_resized), 4);
%     temp = temp_hog(:,:,1:31);

    temp = zeros(scale_model_sz(1), scale_model_sz(2), 1+size(w2c,2));
    if size(im_patch_resized,3) > 1
		gray_patch = rgb2gray(im_patch_resized);
	else
		gray_patch = im_patch_resized;
    end
    % resize it to out size
	gray_patch = imresize(gray_patch, scale_model_sz, 'bilinear', 'AntiAliasing',false);
    % put grayscale channel into output structure
    temp(:, :, 1) = single((gray_patch / 255) - 0.5);
    
    CN = im2c(single(im_patch_resized), w2c, -2);
    CN = imresize(CN, scale_model_sz, 'bilinear', 'AntiAliasing',false);
    % put colornames features into output structure
    temp(:,:,2:(2 + size(w2c, 2) - 1)) = CN;
    
    if s == 1
        out = zeros(numel(temp), nScales, 'single');
    end
    
    % window
    out(:,s) = temp(:) * scale_window(s);
end

end  % endfunction