channelCount=3;
cellCounts=[50 100 150 200];

fid=fopen('simulated_data/sim_cells.txt', 'w');
for cC = 1:length(cellCounts)
    cellCount = cellCounts(cC);
    config = brainbowConfig(channelCount, cellCount);
    cellsToUse=sort(config.cellsToUse);
    fprintf(fid, strcat('simVol_Cell', num2str(cellCount), '_Ch', num2str(channelCount), '.tif\n'));
    for i=1:length(cellsToUse)
        fprintf(fid, num2str(cellsToUse(i)));
        fprintf(fid, '\n');
    end
    [overallRawVolume volumeLabels colorMatrix] = brainbowSimulation_3d_raw(config);

    dim = size(overallRawVolume);
    % 8-bit png
    data = uint8(overallRawVolume*255);
    outputFileName = ['simulated_data/simVol_Cell', num2str(cellCount), '_Ch', num2str(channelCount), '.tif'];
    for i=1:dim(3)
        tmp=squeeze(data(:,:,i,:));
        tmp=permute(tmp, [2, 1, 3]);
        if i==1
            imwrite(tmp, outputFileName);
        else
            imwrite(tmp, outputFileName,'WriteMode','append')
        end
    end
    
    save(['simulated_data/sim_Cell' num2str(cellCount) '_Ch' num2str(channelCount) 'generated_tissue.mat'], 'overallRawVolume','volumeLabels','colorMatrix', '-v7.3');
end

fclose(fid);