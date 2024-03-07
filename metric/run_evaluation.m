% run_evaluation.m
base_path='../';


addpath(fullfile(base_path, 'swc_utils'));
addpath(fullfile(base_path, 'dataset'));

trackers = {'xTracer', 'random'};
result_dirs = {'../dataset/gt', '../dataset/gt'};
gs_dir = '../dataset/gt';
dataset = 'sample1';

metricTypeSet = {'threshold'};
% metricTypeSet = {'trajectory'};

linespecs = config_linespecs;

% Evaluate the results and save to a .mat file
fprintf('\tEvaluating obtained results...\n');
% Create performance mat file for all Thresholds. Later used to plot results
OPE_perfmat(gs_dir, result_dirs, trackers, dataset);

fprintf('\tDrawing Performance plots...\n');
% Rank trackers according to the AUC of the plots
OPE_drawplot(result_dirs{1}, trackers, linespecs, dataset, metricTypeSet);
% TODO: Draw performance plots per challenge of the 

fprintf('\tFinished Evaluation!...\n');