base_path='../';

addpath(fullfile(base_path, 'swc_utils'));

target_swc=load_swc_file('dataset/gt/4-1.swc');
test_swc=load_swc_file('results/exp_real_data/complete/4-1.swc');

target=make_segments(target_swc);
test=make_segments(test_swc);

score=vsa_score(test, target, 10);
rate=BRR(target, test);

fprintf("vsa score: %.2f, BRR: %.2f\n", score, rate)