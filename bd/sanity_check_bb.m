function bb= sanity_check_bb(bb, width, height)
%SANITY_CHECK_BB6 Summary of this function goes here
x=bb(1);
y=bb(2);
w=bb(3);
h=bb(4);

if x < 1
    x=1;
end

if y < 1
    y=1;
end

if x+w > width
    w=width-x;
end

if y+h > height
    y=height-h;
end

if length(bb) == 6
    bb=[x, y, w, h, bb(5), bb(6)];
elseif length(bb) == 4
    bb=[x, y, w, h];
else
    error("check bounding box value, not legal");
end

end

