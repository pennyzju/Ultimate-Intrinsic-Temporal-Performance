function rotated_matrix = rotate_matrix(input_matrix, point1, point2, angle_deg)
    % 计算旋转角度（弧度）
    angle_rad = deg2rad(angle_deg);

    % 计算旋转轴的方向向量并归一化
    axis_vector = point2 - point1;
    axis_vector = axis_vector / norm(axis_vector);

    % 使用 Rodrigues' rotation formula 计算旋转矩阵
    K = [0, -axis_vector(3), axis_vector(2);
         axis_vector(3), 0, -axis_vector(1);
         -axis_vector(2), axis_vector(1), 0];
    rotation_matrix = eye(3) + sin(angle_rad) * K + (1 - cos(angle_rad)) * (K * K);

    % 平移向量，将 point1 平移到原点
    translation_to_origin = -point1;
    translation_back = point1;

    % 获取矩阵的大小
    [m, n, p] = size(input_matrix);

    % 创建新矩阵用于存储旋转后的数据
    rotated_matrix = zeros(m, n, p);

    % 遍历矩阵中的每个点并进行旋转和插值
    for x = 1:m
        for y = 1:n
            for z = 1:p
                % 计算平移后的坐标
                new_coords = [x, y, z] + translation_to_origin;
                % 应用旋转矩阵
                rotated_coords = (rotation_matrix * new_coords')';
                % 平移回去
                final_coords = rotated_coords + translation_back;
                % 插值获取新的值
                new_x = round(final_coords(1));
                new_y = round(final_coords(2));
                new_z = round(final_coords(3));

                % 检查是否在有效范围内
                if new_x > 0 && new_x <= m && new_y > 0 && new_y <= n && new_z > 0 && new_z <= p
                    rotated_matrix(x, y, z) = input_matrix(new_x, new_y, new_z);
                end
            end
        end
    end
end
