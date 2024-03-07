function score_struct=save_eval_xml(score_struct, neuron_id, test_swc_dir, target_swc_dir)
target_swc=load_swc_file(target_swc_dir);
test_swc=load_swc_file(test_swc_dir);

mfd_score=MFD(test_swc(:,3:5), target_swc(:,3:5));

target=make_segments(target_swc);
test=make_segments(test_swc);

score=vsa_score(test, target, 10);
rate=BRR(target, test);

fprintf("Neuron-%s: vsa score: %.2f, MFD: %.2f, BRR: %.2f\n", neuron_id, score, mfd_score, rate)

fn = fieldnames(score_struct);
for i=1:numel(fn)
    if isnumeric(score_struct.(fn{i}).name)
        name=num2str(score_struct.(fn{i}).name);
    else
        name=score_struct.(fn{i}).name;
    end
        
    if( strcmp(name, neuron_id) )
        % update score
        score_struct.(fn{i}).vsa=sprintf('%.4f', score);
        score_struct.(fn{i}).MFD=sprintf('%.2f', mfd_score);
        score_struct.(fn{i}).BRR=sprintf('%.4f', rate);
        score_struct.(fn{i}).length=length(target_swc);
        score_struct.(fn{i}).num_bif=length(target);
    end
end
end