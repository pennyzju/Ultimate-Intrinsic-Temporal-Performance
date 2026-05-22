function [p1, p2] = define_rotation_axis(flag)
    m = 95; n = 112; p = 117;
    switch flag
        case 'left'
            p1 = [0, n/2, 30]; p2 = [m, n/2, 30];
        case 'front'
            p1 = [m/2, 0, 30]; p2 = [m/2, n, 30];
        case 'rot'
            p1 = [m/2, n/2, 0]; p2 = [m/2, n/2, p];
        otherwise
            error('未知的 flag 值');
    end
end