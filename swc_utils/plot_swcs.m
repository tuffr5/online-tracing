function plot_swcs(filenames, legends, titleName, figure_path, figName)
% filenames -cell array

linespecs=config_linespecs();

if ischar(filenames)
    filenames={filenames};
end

h = figure;
p={};
for j=1:size(filenames, 2)
    filename=filenames{1,j};
    % plot swc file
    a=load_swc_file(filename);
    
    child_ids=a(:,1);
    parent_ids=a(:,end);
    ids=setdiff(child_ids, parent_ids);
    ids=[0 ids'];
    
    vx={}; vy={}; vz={};
    for i=1:length(ids)-1
        vx{end+1}=a(ids(i)+1:ids(i+1),3);
        vy{end+1}=a(ids(i)+1:ids(i+1),4);
        vz{end+1}=a(ids(i)+1:ids(i+1),5);
    end
    
    axis tight;
    % use eval to plot multiple lines with same linespec
    base_str_head='plot3(';
    base_str_end= 'linespecs{j})';
    
    exp_str='';
    for i=1:numel(vx)
        exp_str=strcat(exp_str, sprintf('vx{%d}, vy{%d}, vz{%d},', i, i, i));
    end
    p{end+1}=eval(strcat(base_str_head, exp_str,base_str_end));
    hold on;
end

hl=[];
for i=1:numel(p)
    hl= [hl;p{i}(1)];
end

if nargin>1
    legend(hl, legends, 'Location', 'southeast', 'FontSize', 14, 'FontWeight', 'bold');
    box on;
    grid on;
    set(gca,'linewidth',2)
%     title(titleName, 'position', [55,450,110], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
%     set(gca,'xtick',[], 'ytick',[], 'ztick',[]);
%     xlabel(xLabelName);
%     ylabel(yLabelName);
%     saveas(h,fullfile(figure_path, figName), 'png');
    exportgraphics(h, fullfile(figure_path, strcat(figName, '.png')),'Resolution',600); 
    hold off;
end
end


function linespecs = config_linespecs
    linespecs={};
    
    linespecs{end+1}=get_linespec('-', 'r');
    linespecs{end+1}=get_linespec('-', 'b');
    linespecs{end+1}=get_linespec('-', 'm');
    linespecs{end+1}=get_linespec('-', 'y');
    linespecs{end+1}=get_linespec('-', 'c');
    linespecs{end+1}=get_linespec('-', 'g');
    linespecs{end+1}=get_linespec('-', 'k');
end

function linespec = get_linespec(style, color)
    linespec = struct('lineWidth', 3, 'lineStyle', style, 'color', color);
end