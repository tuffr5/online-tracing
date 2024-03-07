function OPE_drawplot(seq_dir, trackers, linespecs, dataset, metricTypeSet)

sequences = dir(fullfile(seq_dir, '*.swc'));
sequences = natsortfiles(sequences);

perfmat_path = fullfile('perfmat', 'OPE');
figure_path = fullfile('figs', 'OPE');

nseq = length(sequences);
ntrk = length(trackers);

% Names of all the trackers evaluated
nameTrkAll = {};
for i = 1:ntrk
    nameTrkAll{end+1} = trackers{i};
end

thresholds = 0:5:100;
rankingType = 'threshold';

% Load the perfmat
perfmat_file = fullfile(perfmat_path, ['perfplot_curves_OPE_' dataset '.mat']);
load(perfmat_file); %'threshold_curve', 'target_length','nameTrkAll'

trajectory_curve={};

for iseq = 1:nseq % for each sequence
    for itrk = 1:ntrk % for each tracker
        trajectory_curve{itrk, iseq} = threshold_curve{itrk, iseq}(1,2);
    end
end

for i=1:length(metricTypeSet)
    metricType = metricTypeSet{i};%'threshold', 'trajectory'
    
    switch metricType
        case 'threshold'
            thresholdSet = thresholds;
            rankIdx = 1;
            curvePlot = threshold_curve;
            figFn = 'threshold_plot';
            titleName = ['Threshold plots of OPE in ' dataset];
            xLabelName = 'Distance threshold';
            yLabelName = 'VSA score';
            
            % cell2array
            curve = reshape(cell2mat(curvePlot), ntrk, length(thresholdSet), nseq);
            curve = squeeze(mean(curve,3));
            
        case 'trajectory'
            thresholdSet = target_length;
            rankIdx = 2;
            curvePlot = trajectory_curve;
            figFn = 'trajectory_plot';
            titleName = ['Trajectory plots of OPE in ' dataset];
            xLabelName = 'Length of gold standard';
            yLabelName = 'VSA score';
            
            % cell2array
            curve = squeeze(reshape(cell2mat(curvePlot), ntrk, nseq));
    end
    
    % Name of image file to save
    figName = [figFn '_OPE_' dataset '_' rankingType];

    % Calculate tracker result for the dataset
    thre = cellfun(@(x)x(rankIdx), curvePlot,'uni',0); 
    perf = mean(cell2mat(thre), 2);


    % Make the legend for the plot with the corresponding ranking type
    for idTrk = 1:ntrk
        trkLegendAll{idTrk} = [nameTrkAll{idTrk} ' [' num2str(perf(idTrk),'%.3f') ']'];
    end
    % Rank the trackers
    [~, trackersRanked] = sort(perf, 'descend');

    h = figure; hold on;
    for idTrk = trackersRanked'
        plot(thresholdSet, curve(idTrk,:), linespecs{idTrk}); hold on;
    end

    legend(trkLegendAll(trackersRanked), 'Location', 'southwest');
    box on; % displays the box outline around the current axes
    title(titleName);
    xlabel(xLabelName);
    ylabel(yLabelName);
    saveas(h,fullfile(figure_path, figName),'png');
        
end

end

