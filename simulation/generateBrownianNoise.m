function voxelColors = generateBrownianNoise(voxelIndices, stackSize, channelCount, opts)

xx_aff                               = zeros(numel(voxelIndices)*13, 1);
yy_aff                               = zeros(numel(voxelIndices)*13, 1);
idx                                  = 1;
for kk = 1:numel(voxelIndices)
  [xx, yy, zz]                       = ind2sub(stackSize, voxelIndices(kk));
  [xxx yyy zzz]                      = meshgrid([-1 0 1], [-1 0 1], [-1 0 1]);
  neighbors                          = repmat([xx yy zz], 27, 1) + [xxx(:) yyy(:) zzz(:)];
  neighbors(14, :)                   = [];
  invalid                            = find(sum(neighbors<1,2) | neighbors(:,1)>stackSize(1) | neighbors(:,2)>stackSize(2) | neighbors(:,3)>stackSize(3));
  neighbors(invalid, :)              = [];
  neighbors                          = sub2ind(stackSize, neighbors(:,1), neighbors(:,2),neighbors(:,3));
  neighbors                          = find(ismember(voxelIndices, neighbors));
  neighbors(neighbors>kk)            = [];
  xx_aff(idx:idx+numel(neighbors)-1) = kk;
  yy_aff(idx:idx+numel(neighbors)-1) = neighbors;
  idx                                = idx+numel(neighbors);
end
xx_aff(idx:end)                      = [];
yy_aff(idx:end)                      = [];
G                                    = sparse(xx_aff, yy_aff, 1, numel(voxelIndices), numel(voxelIndices));
[S, C]                               = graphconncomp(G, 'Weak', true);
representatives                      = zeros(1, S);
disc                                 = [];
for kk = 1:S
  representatives(kk)                = find(C==kk, 1);
  [thisDisc, ~, ~]                   = graphtraverse(G, representatives(kk) , 'Directed', false, 'Method', 'BFS');
  disc                               = [disc thisDisc];
end

voxelColors                          = zeros(numel(voxelIndices), channelCount);
preassigned                          = union(representatives, randi(numel(voxelIndices), 1, round(opts.preAssignedRatio*numel(voxelIndices))));
voxelColors(preassigned, :)          = zeros(numel(preassigned), channelCount);
processed                            = false(size(voxelIndices));
processed(preassigned)               = true;
for kk = 2:numel(voxelIndices)
  if ~ismember(disc(kk), preassigned)

    [xx, yy, zz]                       = ind2sub(stackSize, voxelIndices(disc(kk)));
    [xxx yyy zzz]                      = meshgrid([-1 0 1], [-1 0 1], [-1 0 1]);
    neighbors                          = repmat([xx yy zz], 27, 1) + [xxx(:) yyy(:) zzz(:)];
    neighbors(14, :)                   = [];
    invalid                            = find(sum(neighbors<1,2) | neighbors(:,1)>stackSize(1) | neighbors(:,2)>stackSize(2) | neighbors(:,3)>stackSize(3));
    neighbors(invalid, :)              = [];
    neighbors                          = sub2ind(stackSize, neighbors(:,1), neighbors(:,2), neighbors(:,3));
    neighbors                          = find(ismember(voxelIndices, neighbors));
    neighbors(~processed(neighbors))   = [];

    voxelColors(disc(kk), :)           = mean(voxelColors(neighbors, :), 1) + opts.randomWalkSD*randn(1, channelCount);
    processed(disc(kk))                = true;
  end
end
