function OPE_perfmat(gs_dir, result_dirs, trackers, dataset)

sequences = dir(fullfile(result_dirs{1}, '*.swc'));
sequences = natsortfiles(sequences);

perfmat_path = fullfile('perfmat', 'OPE');

nseq = length(sequences);
ntrk = length(trackers);

% Names of all the trackers evaluated
nameTrkAll = {};
for i = 1:ntrk
    nameTrkAll{end+1} = trackers{i};
end

threshold_curve = cell(ntrk, nseq);
target_length = cell(ntrk, nseq);

for iseq = 1:nseq % for each sequence
    for itrk = 1:ntrk % for each tracker
        fprintf('%-12s,%3d_%-20s\n', sequences(iseq).name,itrk,trackers{itrk});
        
        [scores, target_len] = perf(sequences(iseq), result_dirs{itrk}, gs_dir);
        threshold_curve{itrk, iseq} = scores;
        target_length{itrk, iseq} = target_len;
    end
end

target_length=cell2mat(target_length(1,:));

% save result file
perfmat_file = fullfile(perfmat_path, ['perfplot_curves_OPE_' dataset '.mat']);
save(perfmat_file,'threshold_curve','target_length','nameTrkAll');

end % function OPE_perfmat

function [scores, target_len] = perf(s, results_path, gt_dir)
    % Threshold sampled in the plots
    thresholds = 0:5:100;

    % Load results
    results_file = fullfile(results_path, s.name);
    test_swc=load_swc_file(results_file);
    target_swc=load_swc_file(fullfile(gt_dir, s.name));
    
    target=make_segments(target_swc);
    test=make_segments(test_swc);
    
    target_len=length(target_swc);

    scores=vsa_score(test, target, thresholds);
end % function perf

