function [B, count] = shift_gtone_region(A, dir, d)
    % 找出大于1的索引
    idx = find(A > 1);
    count = numel(idx);  % 大于1的元素个数
    fprintf('大于1的数据个数：%d\n', count);

    % 获取对应的坐标
    [x, y, z] = ind2sub(size(A), idx);

    % 平移坐标
    switch dir
        case 'X'
            x_new = x + d; y_new = y; z_new = z;
        case 'Y'
            x_new = x; y_new = y + d; z_new = z;
        case 'Z'
            x_new = x; y_new = y; z_new = z + d;
        otherwise
            error('方向必须是 ''X'', ''Y'', 或 ''Z''');
    end

    % 去掉越界的点
    valid = x_new >= 1 & x_new <= size(A,1) & ...
            y_new >= 1 & y_new <= size(A,2) & ...
            z_new >= 1 & z_new <= size(A,3);

    x_new = x_new(valid);
    y_new = y_new(valid);
    z_new = z_new(valid);
    values = A(sub2ind(size(A), x(valid), y(valid), z(valid)));

    % 构造新矩阵
    B = zeros(size(A));
    inds_new = sub2ind(size(A), x_new, y_new, z_new);
    B(inds_new) = values;
end
