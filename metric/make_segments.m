function segments=make_segments(a)
segments={};

% find end_points with no children
child_ids=a(:,1);
parent_ids=a(:,end);
ids=setdiff(child_ids, parent_ids);

ids=[0 ids'];

ls=max(size(ids));

if ls==1
    segments{1,1}=a(:,3:5);
else
    for i=1:ls-1
        try
            segments{i,1}=a(ids(i)+1:ids(i+1),3:5);
        catch
            error("Error.")
        end
    end
    
    if ids(end)<size(a,1)
        segments{i+1,1}=a(ids(i+1):end,3:5);
    end
end
end