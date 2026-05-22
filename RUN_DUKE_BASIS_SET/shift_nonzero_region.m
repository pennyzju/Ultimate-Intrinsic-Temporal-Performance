function B = shift_nonzero_region(A, dir, d)
    % 找到数据的有效范围（非零区域）
    [x, y, z] = ind2sub(size(A), find(A));  
    xmin = min(x); xmax = max(x);
    ymin = min(y); ymax = max(y);
    zmin = min(z); zmax = max(z);

    % 提取非零区域数据
    A_crop = A(xmin:xmax, ymin:ymax, zmin:zmax);
    [Xsize, Ysize, Zsize] = size(A_crop);

    % 确定平移方向
    new_xmin = xmin; new_ymin = ymin; new_zmin = zmin;
    switch dir
        case 'X', new_xmin = xmin + d;
        case 'Y', new_ymin = ymin + d;
        case 'Z', new_zmin = zmin + d;
        otherwise, error('方向必须是 "X", "Y", 或 "Z"');
    end

    % 创建与 A 相同大小的零数组
    B = zeros(size(A));

    % 计算新数据的范围，确保不超出边界
    x_range = max(1, new_xmin) : min(size(A, 1), new_xmin + Xsize - 1);
    y_range = max(1, new_ymin) : min(size(A, 2), new_ymin + Ysize - 1);
    z_range = max(1, new_zmin) : min(size(A, 3), new_zmin + Zsize - 1);

    % 计算裁剪数据的有效索引（防止超出范围）
    x_crop = 1:length(x_range);
    y_crop = 1:length(y_range);
    z_crop = 1:length(z_range);

    % 填充 B
    B(x_range, y_range, z_range) = A_crop(x_crop, y_crop, z_crop);
end

