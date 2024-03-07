function score = vsa_score(test, target, thresholds)
%calculation for variation of segments alignment score
% both inputs should be the form of Nx1 cell array

for i=1:size(target, 1)
    if size(target{i},1)<=10
        target{i}=[];
    end
end

% for i=1:size(test, 1)
%     if size(test{i},1)<=10
%         test{i}=[];
%     end
% end

target=target(~cellfun('isempty',target));
% test=test(~cellfun('isempty',test));

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


match_l=zeros(lT, lS);
for i=1:lT
    for j=1:lS
        match_l(i,j)=count_matches(test{i}, target{j}, 15);
    end
end
%% step 2: local alignment for each pair and count matches
% one test segment can map to multiple targets, but not the reverse
score_m=zeros(lS, length(thresholds));

for i=1:lS
    tmp=match_l(:,i);
    j=find(tmp==max(tmp));
    j=j(1);
    test_segs_len=size(test{j}, 1);
    matches=count_matches(test{j}, target{i}, thresholds);
    score_m(i,:)=min(matches/seg_len(i), matches/size(test{j},1));
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