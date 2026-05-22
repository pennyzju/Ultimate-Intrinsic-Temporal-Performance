flag = 'front';
angle_deg = -5;
% 获取数据维度
    m = 95;
    n = 112;
    p = 117;

    % 根据flag的值选择旋转轴
    if strcmp(flag, 'left')
        point1 = [0, n/2, 30];
        point2 = [m, n/2, 30];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    elseif strcmp(flag, 'front')
        point1 = [m/2, 0, 30];
        point2 = [m/2, n, 30];
    elseif strcmp(flag, 'rot')
        point1 = [m/2, n/2, 0];
        point2 = [m/2, n/2, p];
    else
        error('Invalid flag value. Flag should be ''left'', ''front'', or ''rot''.');
    end

    mask2 = rotate_matrix(mask, point1, point2, angle_deg);