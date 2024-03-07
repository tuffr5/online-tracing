clear all;
clc;

base_dir='/path/tp/xBTracer/';
other_trackers_dir='/path/tp/xBTracer/baseline_trackers/';

fprintf('Part 1: Run comparison experiments\n');

cd(base_dir);

trackers={'xBTracer', 'csr-dcf', 'STRCF', 'ASRCF', 'AutoTrack'};

for i=1:numel(trackers)
	if ~strcmp('xBTracer', trackers{i})
		cd(strcat(other_trackers_dir, trackers{i}, '/'));
	end
	fprintf('Run %s on real-world Dataset:\n', trackers{i});
	exp_name=strcat(trackers{i}, '_real_data/');
	run_on_realData(exp_name);
	fprintf('Finish %s on real-world Dataset.\n', trackers{i});
	fprintf('Run %s on simulated Dataset:\n', trackers{i});
	exp_name=strcat(trackers{i}, '_sim_data/');
	run_on_simData(exp_name);
	fprintf('Finish %s on simulated Dataset.\n', trackers{i});
end


fprintf('Part 2: Run ablation experiments\n');

cd(base_dir);

ablations={'wo_c', 'wo_b', 'wo_m'};

for i=1:numel(ablations)
	if strcmp('wo_c', ablations{i})
		config.enable_bifurcation=true;
		config.enable_cross_plane=false;
		config.enable_alpha_mask=true;
	end

	if strcmp('wo_b', ablations{i})
		config.enable_bifurcation=false;
		config.enable_cross_plane=true;
		config.enable_alpha_mask=true;
	end 

	if strcmp('wo_m', ablations{i})
		config.enable_bifurcation=true;
		config.enable_cross_plane=true;
		config.enable_alpha_mask=false;
	end 
	
	fprintf('Run %s on real-world Dataset:\n', ablations{i});
	exp_name=strcat('xBTracer_real_data_', ablations{i}, '/');
	run_on_realData(exp_name, config);
	fprintf('Finish %s on real-world Dataset.\n', ablations{i});
	fprintf('Run %s on simulated Dataset:\n', ablations{i});
	exp_name=strcat('xBTracer_sim_data_', ablations{i}, '/');
	run_on_simData(exp_name, config);
	fprintf('Finish %s on simulated Dataset.\n', ablations{i});
end

fprintf('Congratulations!!! Finish all experiments!!!\n');