function regions=get_bb(X, Y)  
    x1 = round(min(X));
    x2 = round(max(X));
    y1 = round(min(Y));
    y2 = round(max(Y));
    regions = round([x1, y1, x2 - x1 + 1, y2 - y1 + 1]);
end