tif_dir='simulation/simulated_data/*.tif';

tifs=dir(tif_dir);
tifs=natsortfiles(tifs);

for i=1:numel(tifs)
    tif=fullfile(tifs(i).folder, tifs(i).name);
    protein_density(tif);
end

real='/home/path/to/Downloads/nTracer_sample_img/nTracer sample.tif';
protein_density(real);


function protein_density(tif, threshold)
if nargin < 2
    threshold=50;
end

[imgs,h,w]=imread_big(tif);
%Height*Width*(CZ) --> Height*Width*C*Z !!!important
imgs = reshape(imgs,h,w,3,[]);
% YXCZ -> ZYXC
imgs=permute(imgs,[4 1 2 3]);
imgs_sum=sum(imgs, 4);
mask_imgs=imgs_sum>threshold;

ratio=sum(mask_imgs(:))/numel(imgs_sum)*100;
disp(ratio);

end

    