function rate = BRR(target, test)
% Bifurcation retrival rate
for i=1:size(target, 1)
    if size(target{i},1)<=20
        target{i}=[];
    end
end

% for i=1:size(test, 1)
%     if size(test{i},1)<=10
%         test{i}=[];
%     end
% end

target=target(~cellfun('isempty',target));
test=test(~cellfun('isempty',test));
if size(target,1) <= 1
    rate=0;
else
    rate=min((size(test,1)-1)/(size(target,1)-1), 1);
end
end

