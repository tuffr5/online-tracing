function xyz=wrap_xyz(c, plane, idx)
    if plane==1
        xyz=[c(1) c(2) idx];
    end

    if plane==2
        xyz=[idx c(2) c(1)];
    end

    if plane==3
        xyz=[c(1) idx c(2)];
    end 
end