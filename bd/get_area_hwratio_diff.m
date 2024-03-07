function [diff_area, diff_hw_ratio] = get_area_hwratio_diff(bb, base_area, base_hwratio)
%bb 4X2
h=norm(bb(1,:)-bb(2,:));
w=norm(bb(2,:)-bb(3,:));
    
area= h*w;
hw_ratio=max(h,w)/min(h,w);

if base_area ~= 0
    diff_area=abs(area-base_area)/base_area;
else
    diff_area = area;
end

if base_hwratio ~= 0
    diff_hw_ratio=abs(hw_ratio-base_hwratio)/base_hwratio;
else
    diff_hw_ratio=hw_ratio;
end
end

