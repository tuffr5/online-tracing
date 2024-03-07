function [path, cost]=Astar(MAP, start, goal)
%Modified from the code of Anthony Chrabieh

%Heuristic Weight%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight = 2.5; %Try 1, 1.5, 2, 2.5
%Increasing weight makes the algorithm greedier, and likely to take a
%longer path, but with less computations.
%weight = 0 gives Djikstra algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=size(MAP);
H=inf*ones(a(1:2));
G=inf*ones(a(1:2));

mean_intens=mean_of_wins([start; goal], MAP);
start_color=mean_intens(1,:);% 3x1
goal_color=mean_intens(2,:);% 3x1

%Heuristic Map of all nodes
for x = 1:size(MAP,1)
    for y = 1:size(MAP,2)
        if(MAP(x,y)~=inf)
            xy_color=mean_of_wins([x, y], MAP);
            color_dist_1=abs(start_color-xy_color)./start_color;
            color_dist_2=abs(goal_color-xy_color)./goal_color;
            
            H(x,y) = sum(abs(goal-[x,y]), 'all')...
                + weight*min(sum(color_dist_1(:)), sum(color_dist_2(:))); % manhattan distance
        end
    end
end

%initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G(start(1),start(2)) = 0;
F(start(1),start(2)) = H(start(1),start(2));

closedNodes = [];
openNodes = [start G(start(1),start(2)) F(start(1),start(2)) 0]; %[x y G F cameFrom]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve
solved = false;
while(~isempty(openNodes))
    %find node from open set with smallest F value
    [A,I] = min(openNodes(:,4));
    
    %set current node
    current = openNodes(I,:);
    
    %if goal is reached, break the loop
    if(current(1:2)==goal)
        closedNodes = [closedNodes;current];
        solved = true;
        break;
    end
    
    %remove current node from open set and add it to closed set
    openNodes(I,:) = [];
    closedNodes = [closedNodes;current];
    
    %for all neighbors of current node
    for x = current(1)-1:current(1)+1
        for y = current(2)-1:current(2)+1
            
            %if out of range skip
            if (x<1||x>size(MAP,1)||y<1||y>size(MAP,2))
                continue
            end
            
            %if current node skip
            if (x==current(1)&&y==current(2))
                continue
            end
            
            %if already in closed set skip
            skip = 0;
            for j = 1:size(closedNodes,1)
                if(x == closedNodes(j,1) && y==closedNodes(j,2))
                    skip = 1;
                    break;
                end
            end
            if(skip == 1)
                continue
            end
            
            A = [];
            %Check if already in open set
            if(~isempty(openNodes))
                for j = 1:size(openNodes,1)
                    if(x == openNodes(j,1) && y==openNodes(j,2))
                        A = j;
                        break;
                    end
                end
            end
            
            newG = G(current(1),current(2)) ...
                + sum(abs([current(1) ,current(2)]-[x, y]),'all');% manhattan distance
            
            %if not in open set, add to open set
            if(isempty(A))
                G(x,y) = newG;
                newF = G(x,y) + H(x,y);
                newNode = [x y G(x,y) newF size(closedNodes,1)];
                openNodes = [openNodes; newNode];
            end
            
            %if no better path, skip
            if (newG >= G(x,y))
                continue
            end
            
            G(x,y) = newG;
            newF = newG + H(x,y);
            openNodes(A,3:5) = [newG newF size(closedNodes,1)];
        end
    end
end

if (solved)
    %Path plotting
    j = size(closedNodes,1);
    path = [];
    low_intensitys=0;
    up_intensitys=0;
    while(j > 0)
        x = closedNodes(j,1);
        y = closedNodes(j,2);
        j = closedNodes(j,5);
        path = [x,y;path];
        
        xy_color=mean_of_wins([x, y], MAP);
        color_dist_1=abs(start_color(1)-xy_color(1))./start_color(1); % only use hue
        color_dist_2=abs(goal_color(1)-xy_color(1))./goal_color(1);
        
        up_intensitys=up_intensitys+max(sum(color_dist_1(:)), sum(color_dist_2(:)));
        low_intensitys=low_intensitys+min(sum(color_dist_1(:)), sum(color_dist_2(:)));
    end
    cost=0.5*(up_intensitys+low_intensitys)/length(path);
    
else
    path = [];
    cost = inf;
end
end

function mean_intens=mean_of_wins(vals, MAP, win_sz)
if nargin == 2
    win_sz = 2;
end

mean_intens=[];
for i =1:size(vals, 1)
    mask_win=zeros(size(MAP, 1), size(MAP, 2));
    mask_win(vals(i,1), vals(i,2))=1;
    % get a window by dilation
    mask_win=imdilate(mask_win, strel('square', win_sz));

    img_mask=MAP .* double(mask_win);
    mean_inten=sum(img_mask, [1 2])/sum(mask_win(:));
    mean_intens=[mean_intens; squeeze(mean_inten)'];
end

end