% SWC format: http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
function save_to_swc_file(tt, outfilename)
if isempty(tt)
    return;
end

if size(tt,2)<7
    error('The first variable must have at least 7 columns.'),
end

f = fopen(outfilename, 'w');
if f<=0
    error('Fail to open file to write');
end

for i=1:size(tt,1)
    fprintf(f, '%d %d %d %d %d %5.3f %d\n', tt(i,1), tt(i,2), tt(i,3), tt(i,4), tt(i,5), tt(i,6), tt(i,7));
end

fclose(f);
end