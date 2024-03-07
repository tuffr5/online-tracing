function [overallRawVolume volumeLabels colorMatrix] = brainbowSimulation_3d_raw(config)
%function [overallRawVolume volumeLabels colorMatrix noiseImage] = brainbowSimulation_3d(config)

cellsToUse              = config.cellsToUse;
ANGLESTEP               = config.ANGLESTEP;
SHIFTSTEP               = config.SHIFTSTEP;
ZSHIFTSTEP              = config.ZSHIFTSTEP;
brightnessThreshold     = 0.5;
normalSD                = config.normalSD;
xSize                   = config.xSize;
ySize                   = config.ySize;
zSize                   = config.zSize;
overflow                = config.overflow;
channelCount            = config.channelCount;
overallRawVolume        = zeros(xSize, ySize, zSize, channelCount);
overallTraceBasedVolume = zeros(xSize, ySize, zSize);
stackSize               = [xSize ySize zSize];
code = cellListArray; code = code(cellsToUse);
cutLocs = cutLocationsArray; cutLocs = cutLocs(cellsToUse,:);
allSizes = zeros(numel(code),2);
colorMatrix = pickRandomColors(numel(cellsToUse), channelCount);
% INSERT THE FIRST CELL AS IS
kk      = 1;
retNo   = num2str(code{kk}{1}); cellNo = num2str(code{kk}{2}); cellID = ['Image',retNo,cellNo,'_01'];

info    = imfinfo(['simulated_data/material/',num2str(code{kk}{1}),'-',num2str(code{kk}{2}),'.tif']);

thisVol = zeros(info(1).Height, info(1).Width, numel(info));
for zz = 1:numel(info);
  thisVol(:,:,zz) = imread(['simulated_data/material/',num2str(code{kk}{1}),'-',num2str(code{kk}{2}),'.tif'],'Index',zz);
end;
thisVol = (thisVol(cutLocs(kk,1):cutLocs(kk,2), cutLocs(kk,3):cutLocs(kk,4), cutLocs(kk,5):cutLocs(kk,6))>0); [xSize_tV ySize_tV zSize_tV] = size(thisVol);
cuts = [max(1, round((xSize_tV-xSize)/2+1)) min(xSize_tV, round((xSize_tV+xSize)/2)) max(1, round((ySize_tV-ySize)/2+1)) min(ySize_tV, round((ySize_tV+ySize)/2)) max(1, round((zSize_tV-zSize)/2+1)) min(zSize_tV, round((zSize_tV+zSize)/2))];
thisVol = thisVol(cuts(1):cuts(2), cuts(3):cuts(4), cuts(5):cuts(6));
%thisConnComp = bwconncomp(thisVol); toRemove = []; for oo=config.maxConnCompPerCell+1:thisConnComp.NumObjects; toRemove = [toRemove; thisConnComp.PixelIdxList{oo}]; end; thisVol(toRemove) = false;
overallTraceBasedVolume(1:size(thisVol,1), 1:size(thisVol,2), 1:size(thisVol,3)) = thisVol;
volumeLabels{1} = false(xSize, ySize, zSize);
volumeLabels{1}(1:size(thisVol,1), 1:size(thisVol,2), 1:size(thisVol,3)) = thisVol;
overallCounter  = overallTraceBasedVolume;
voxelIndices = find(volumeLabels{1});
voxelColors = generateBrownianNoise(voxelIndices, stackSize, channelCount, config.colors);
for ch = 1:channelCount; overallRawVolume(voxelIndices+prod(stackSize)*(ch-1)) = colorMatrix(1,ch) + voxelColors(:,ch); end;
% TRY TO MINIMIZE THE OVERLAP IN PLACING THE REMAINING CELLS
for kk = 2 : numel(code)
  retNo   = num2str(code{kk}{1}); cellNo = num2str(code{kk}{2}); cellID = ['Image',retNo,cellNo,'_01'];
  info    = imfinfo(['simulated_data/material/',num2str(code{kk}{1}),'-',num2str(code{kk}{2}),'.tif']);
  thisVol = zeros(info(1).Height, info(1).Width, numel(info));
  for zz = 1:numel(info);
    thisVol(:,:,zz) = imread(['simulated_data/material/',num2str(code{kk}{1}),'-',num2str(code{kk}{2}),'.tif'],'Index',zz);
  end;
  thisVol = (thisVol(cutLocs(kk,1):cutLocs(kk,2), cutLocs(kk,3):cutLocs(kk,4), cutLocs(kk,5):cutLocs(kk,6))>0); [xSize_tV ySize_tV zSize_tV] = size(thisVol);
  cuts = [max(1, round((xSize_tV-xSize)/2+1)) min(xSize_tV, round((xSize_tV+xSize)/2)) max(1, round((ySize_tV-ySize)/2+1)) min(ySize_tV, round((ySize_tV+ySize)/2)) max(1, round((zSize_tV-zSize)/2+1)) min(zSize_tV, round((zSize_tV+zSize)/2))];
  thisVol= thisVol(cuts(1):cuts(2), cuts(3):cuts(4), cuts(5):cuts(6));
  thisConnComp = bwconncomp(thisVol); toRemove = []; for oo=config.maxConnCompPerCell+1:thisConnComp.NumObjects; toRemove = [toRemove; thisConnComp.PixelIdxList{oo}]; end;thisVol(toRemove) = false;
  minOverlap = prod(size(overallRawVolume));
  for rotAngle = 0 : ANGLESTEP : 360-ANGLESTEP
    vol = imrotate(thisVol,rotAngle,'nearest','crop');
    [thisXsize, thisYsize, thisZsize] = size(vol);
    % SHIFT VALUES ARE SORTED TO EVALUATE THE OVERLAPS OF NO-OVERFLOW CASES LAST
    allXshifts = [-overflow : SHIFTSTEP : xSize-thisXsize+overflow];
    allXshifts = [allXshifts(allXshifts<0) allXshifts(allXshifts>xSize-thisXsize) allXshifts(allXshifts>=0 & allXshifts<=xSize-thisXsize)];
    allYshifts = [-overflow : SHIFTSTEP : ySize-thisYsize+overflow];
    allYshifts = [allYshifts(allYshifts<0) allYshifts(allYshifts>ySize-thisYsize) allYshifts(allYshifts>=0 & allYshifts<=ySize-thisYsize)];
    % Z JITTER IS SMALLER FOR MORE REALISTIC ASSEMBLY OF RGCs AND TO INCREASE THE NUMBER OF OVERLAPS
    allZshifts = [max(-5, -overflow) : ZSHIFTSTEP : min(5, overflow)]; %allZshifts = [-overflow : SHIFTSTEP : zSize-thisZsize+overflow];
    allZshifts = [allZshifts(allZshifts<0) allZshifts(allZshifts>zSize-thisZsize) allZshifts(allZshifts>=0 & allZshifts<=zSize-thisZsize)];
    for xcounter = 1 : numel(allXshifts)
      xShift = allXshifts(xcounter);
      for ycounter = 1 : numel(allYshifts)
        yShift = allYshifts(ycounter);
        for zcounter = 1 : numel(allZshifts)
          zShift = allZshifts(zcounter);
          tmpVol = false(xSize, ySize, zSize);
          xl = max(1, 1+xShift); xh = min(xSize, thisXsize+xShift);
          yl = max(1, 1+yShift); yh = min(ySize, thisYsize+yShift);
          zl = max(1, 1+zShift); zh = min(zSize, thisZsize+zShift);
          tmpVol(xl:xh, yl:yh, zl:zh) = vol(xl-xShift:xh-xShift, yl-yShift:yh-yShift, zl-zShift:zh-zShift);
          thisOverlap = nnz(tmpVol & overallTraceBasedVolume);
          if thisOverlap <= minOverlap
            minOverlap = thisOverlap;
            bestTraceBasedVolume = tmpVol;
          end
        end
      end
    end
  end
  overallTraceBasedVolume = overallTraceBasedVolume | bestTraceBasedVolume;
  volumeLabels{kk}        = bestTraceBasedVolume;
  overallCounter          = overallCounter + bestTraceBasedVolume;
  voxelIndices            = find(bestTraceBasedVolume);
  voxelColors             = generateBrownianNoise(voxelIndices, stackSize, channelCount, config.colors);
  for ch = 1:channelCount; overallRawVolume(voxelIndices+prod(stackSize)*(ch-1)) = colorMatrix(kk,ch) + voxelColors(:,ch); end;
end
% AVERAGE WHERE NEURITES COEXIST
overallRawVolume(:,:,:,1) = overallRawVolume(:,:,:,1)./overallCounter;
overallRawVolume(:,:,:,2) = overallRawVolume(:,:,:,2)./overallCounter;
overallRawVolume(:,:,:,3) = overallRawVolume(:,:,:,3)./overallCounter;
overallRawVolume(isnan(overallRawVolume)) = 0;
% ADD POISSON NOISE
overallRawVolume(:,:,:,1) = overallRawVolume(:,:,:,1) + randn(xSize, ySize, zSize) * normalSD;
overallRawVolume(:,:,:,2) = overallRawVolume(:,:,:,2) + randn(xSize, ySize, zSize) * normalSD;
overallRawVolume(:,:,:,3) = overallRawVolume(:,:,:,3) + randn(xSize, ySize, zSize) * normalSD;
% MAXIMUM VALUE IN EACH READOUT COLOR IS 1 - OTHERWISE COLORS SATURATE AT 1
overallRawVolume = max( min( overallRawVolume, 1 ), 0 );
%noiseImage       = overallRawVolume - noiseFreeVolume;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colorMatrix = pickRandomColors(count, chCount, brightnessThreshold)

if nargin<3; brightnessThreshold = 0; end;

colorMatrix = zeros(count, chCount);
for neuron = 1:count
  randColor = zeros(1, chCount);
  for ch=1:chCount; randColor(ch) = rand; end;
  while max(randColor)<brightnessThreshold
    for ch=1:chCount; randColor(ch) = rand; end;
  end
  colorMatrix(neuron, :) = randColor;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code = cellListArray
code = cell(0);
materials=dir(['simulated_data/material/*.tif']);
for i=1:numel(materials)
  [retNo, cellNo]=parse_name(materials(i).name);
  code{end+1}{1} = retNo; code{end}{2} = cellNo;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutLocs = cutLocationsArray
cutLocs=[];

materials=dir(['simulated_data/material/*.tif']);
for i=1:numel(materials)
    info=imfinfo(['simulated_data/material/',materials(i).name]);
    cutLocs=[cutLocs; [1, info(1).Height,1, info(1).Width, 1, numel(info)]];
end
end

function [retNo, cellNo]=parse_name(name)
str_splited=split(name, '-');
retNo=str_splited(1);
retNo=retNo{1};

remaining=str_splited(2);
remaining_splited=split(remaining, '.tif');

cellNo=remaining_splited(1);
cellNo=cellNo{1};
end