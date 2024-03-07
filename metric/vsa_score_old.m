function score = vsa_score(test, target, thresholds)
%calculation for variation of segments alignment score
% both inputs should be the form of Nx1 cell array

if nargin==2
    thresholds = 10;
end

%% step 1: identify pairs of test and target segments
lT=size(test, 1);
lS=size(target, 1);

seg_len=[];
target_len=0;

for i=1:lS
    seg_len=[seg_len size(target{i},1)];
    target_len=target_len+size(target{i},1);
end


dist_l=zeros(lT, lS);
for i=1:lT
    for j=1:lS
        dist_l(i,j)=DiscreteFrechetDist(test{i}, target{j});
    end
end

[C, col]=min(dist_l, [], 2);

%% step 2: local alignment for each pair and count matches
% one test segment can only map to one target
score_m=zeros(lS, length(thresholds));
cols=unique(col);

for i=1:max(size(cols))
    rows=find(col==cols(i));
    
    test_segs_len=0;
    tmp=[];
    
    for j=1:max(size(rows))
        
        if isinf(C(rows(j)))
           tmp=[];
        else
           tmp=[tmp; test{rows(j)}];
        end   
    end
    
    if ~isempty(tmp)
        matches=count_matches(tmp, target{cols(i)}, thresholds);
        test_segs_len=test_segs_len+size(tmp, 1);

        score_m(cols(i),:)=0.9*matches/seg_len(cols(i)) + 0.1*matches/test_segs_len;
    end
    
end
%% step 3: calculate the final score proportinal to the len of total target
score=zeros(1,length(thresholds));
for i=1:lS
    score=score+seg_len(i)/target_len*score_m(i,:);
end

end

function matches=count_matches(s1, s2, thresholds)
%% count matches by thresholding
dist=pdist2(s1, s2, 'euclidean');
matches=zeros(length(thresholds),1);

for i=1:length(thresholds)
    dist=dist<=thresholds(i)*sqrt(3);

    % make sure one data points in target tracing maximumly mapped one time
    % using Graph Matching (rank of the dist matrix)
    A=sparse(dist);
    matches(i)=sprank(A);
end

end