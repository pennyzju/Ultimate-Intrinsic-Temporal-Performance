function [point1, point2] = get_rotation_axis(flag)
    % 获取旋转轴的两个端点
    % 输入:
    %   flag: 'left', 'front' 或 'rot'
    % 输出:
    %   point1, point2: 两个三维点，定义旋转轴

    % 数据维度
    m = 95;
    n = 112;
    p = 117;

    % 根据flag的值选择旋转轴
    switch flag
        case 'left'
            point1 = [0, n/2, 30];
            point2 = [m, n/2, 30];

        case 'front'
            point1 = [m/2, 0, 30];
            point2 = [m/2, n, 30];

        case 'rot'
            point1 = [m/2, n/2, 0];
            point2 = [m/2, n/2, p];

        otherwise
            error('Invalid flag value. Flag should be ''left'', ''front'', or ''rot''.');
    end
end
