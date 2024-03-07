function dist = MFD(test, target)

distance=pdist2(test, target);
min_distance=min(distance, [], 2);

dist=mean(min_distance(:));
end

