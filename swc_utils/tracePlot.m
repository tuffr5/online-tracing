function tracePlot(data)
%--------------------------------------------------------------------------
% Date: July 4, 2020
% input: tracing result
% output:
%---------------------------------------------------------------------------
neurons = unique([data{:,1}]);
count = length(neurons);

color = hsv(count);

figure;

for j = 1:count
    idx = find(cell2mat(data(:,1)) == neurons(j));
    for i = idx(1):idx(length(idx))
        x = data{i,3}(:,1);
        y = data{i,3}(:,2);
        z = data{i,3}(:,3);
%         set(gca, 'ydir', 'reverse', 'zdir', 'reverse', 'Color', 'k', 'xtick', [], 'ytick', [], 'ztick', [])
        plot3(x, y, z, 'Color', color(j,:))
        hold on;
    end
end
end
