function plot_swcs(filenames)
% filenames -cell array

color = 'brgmcyk';

figure;
tiledlayout(size(filenames, 2),1);

if ischar(filenames)
    filenames={filenames};
end

for j=1:size(filenames, 2)
    filename=filenames{1,j};
    % plot swc file
    a=load_swc_file(filename);
    
    ids= find(a(:,1) ~= a(:,end)+1);
    
    if ~isempty(ids)
        if ids(1)==1
            ids(1)=[];
        end
    end

    
    if ~isempty(ids) && sum(ids~=1)==size(ids,1)
        if size(ids, 1) > 1
            for i=1:size(ids,1)-1
                plot3(a([a(ids(i),end), ids(i):ids(i+1)-1],3), ...
                    a([a(ids(i),end), ids(i):ids(i+1)-1],4), ...
                    a([a(ids(i),end), ids(i):ids(i+1)-1],5), ...
                    color(j));
                hold on;
            end
        end

        plot3(a(1:ids(1)-1,3), a(1:ids(1)-1,4), a(1:ids(1)-1,5), color(j));
        hold on;
        plot3(a([a(ids(end),end),ids(end):end],3),...
            a([a(ids(end),end),ids(end):end],4),...
            a([a(ids(end),end),ids(end):end],5),...
            color(j));  
        hold on;
    else
        plot3(a(:,3), a(:,4), a(:,5), color(j));
        hold on;
    end
end

hold off;
end