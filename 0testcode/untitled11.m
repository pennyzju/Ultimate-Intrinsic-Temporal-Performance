% 假设三维数组 data，维度为 (X, Y, Z)，理论上 Z = Z_mid 存在数据
data = usnr_map_anti;
% 假设三维数组 data，维度为 (X, Y, Z)
[X, Y, Z] = size(data);

% 计算三个方向的中间层
X_mid = ceil(X / 2);
Y_mid = ceil(Y / 2);
Z_mid = ceil(Z / 2);

% 计算三个中间层的数据特征（众数或均值）
mid_X_layer = squeeze(data(X_mid, :, :));
mid_Y_layer = squeeze(data(:, Y_mid, :));
mid_Z_layer = squeeze(data(:, :, Z_mid));

% 去除零值（假设零值为无效数据）
valid_X = mid_X_layer(mid_X_layer ~= 0);
valid_Y = mid_Y_layer(mid_Y_layer ~= 0);
valid_Z = mid_Z_layer(mid_Z_layer ~= 0);

% 计算参考值（使用众数 mode 或均值 mean）
if ~isempty(valid_X), ref_X = mode(valid_X(:)); else, ref_X = 0; end
if ~isempty(valid_Y), ref_Y = mode(valid_Y(:)); else, ref_Y = 0; end
if ~isempty(valid_Z), ref_Z = mode(valid_Z(:)); else, ref_Z = 0; end

% 遍历所有点，修正错位数据
for x = 1:X
    for y = 1:Y
        for z = 1:Z
            if data(x, y, z) ~= 0 % 非空数据
                % 计算该点到三个中间层的距离
                dx = abs(x - X_mid);
                dy = abs(y - Y_mid);
                dz = abs(z - Z_mid);

                % 判断是否更接近中间层的数据
                moveToX = (dx > 0) && (abs(data(x, y, z) - ref_X) < abs(data(X_mid, y, z) - ref_X));
                moveToY = (dy > 0) && (abs(data(x, y, z) - ref_Y) < abs(data(x, Y_mid, z) - ref_Y));
                moveToZ = (dz > 0) && (abs(data(x, y, z) - ref_Z) < abs(data(x, y, Z_mid) - ref_Z));

                % 依次移动数据到最近的中间层
                if moveToX
                    data(X_mid, y, z) = data(x, y, z);
                    data(x, y, z) = 0;
                elseif moveToY
                    data(x, Y_mid, z) = data(x, y, z);
                    data(x, y, z) = 0;
                elseif moveToZ
                    data(x, y, Z_mid) = data(x, y, z);
                    data(x, y, z) = 0;
                end
            end
        end
    end
end

% 结果数据 data 已修正

visualize_3D_slices_correctly(usnr_map_anti);
visualize_3D_slices_correctly(data);

[X, Y, Z] = meshgrid(1:size(data,2), 1:size(data,1), 1:size(data,3)); 
figure;
slice(X, Y, Z, data, [], [46:1:48], []); % 在 Z=5,10,15 层切片
colormap(jet); 
colorbar;
shading interp;
axis equal;
title('多层切片图');

figure;
threshold_value = 0;
[x, y, z] = ind2sub(size(data), find(data > threshold_value)); % 仅绘制非零点
scatter3(x, y, z, 20, data(sub2ind(size(data), x, y, z)), 'filled');
colormap(jet);
colorbar;
axis equal;
grid on;
title('三维散点图');


% 假设三维数组 data，尺寸为 (X, Y, Z)
[X, Y, Z] = size(data);

% 计算中间层索引
X_mid = ceil(X / 2);
Y_mid = ceil(Y / 2);
Z_mid = ceil(Z / 2);

% 遍历非中间层，查找非零数据
[x_idx, y_idx, z_idx] = ind2sub(size(data), find(data ~= 0)); % 找到所有非零元素索引

% 过滤出位于非中间层的数据点
non_mid_mask = (x_idx ~= X_mid) & (y_idx ~= Y_mid) & (z_idx ~= Z_mid);
x_out = x_idx(non_mid_mask);
y_out = y_idx(non_mid_mask);
z_out = z_idx(non_mid_mask);

% 输出结果
if isempty(x_out)
    fprintf('未发现非中间层上的非零数据。\n');
else
    fprintf('非中间层上的非零数据坐标 (X, Y, Z):\n');
    disp([x_out, y_out, z_out]);
end
data(X_mid, :, :) = 0; % X 方向中间层
data(:, Y_mid, :) = 0; % Y 方向中间层
data(:, :, Z_mid) = 0; % Z 方向中间层

