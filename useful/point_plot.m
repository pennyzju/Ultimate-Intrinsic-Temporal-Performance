%% 1. 数据准备
load('H:\UISNR\20240801_UISNR_output\left_0\default_material_maps.mat')
V = epsilon_r;

% 定义我们要看的切片层数 (即第 1 维度的索引)
target_slice_idx = 48;

% 定义坐标点 [Row(层), Col(行), Page(列)]
% 注意：这里的顺序对应矩阵索引 V(48, 56, 74)
point_list = [
    48, 56, 74;   % 点 1
    48, 60, 74;   % 点 2
    48, 70, 74;   % 点 3
    48, 80, 74    % 点 4
];

%% 2. 提取切片与筛选点
% 提取第 48 层的二维图像
% squeeze 将 (1, 100, 100) 压缩为 (100, 100)
slice_img = squeeze(V(target_slice_idx, :, :));

% 筛选出位于当前层（第48层）的点
% 逻辑：只画那些第一维坐标等于 target_slice_idx 的点
points_on_slice = point_list(point_list(:,1) == target_slice_idx, :);

%% 3. 可视化绘图
figure('Color', 'w');

% (1) 绘制热力图/切片图
imagesc(slice_img);
colormap('gray'); % MRI 通常用 gray 或 jet
colorbar;
axis image;       % 保持长宽比一致，防止图像变形

hold on; % 保持图像，准备叠加点

% (2) 叠加标记点
if ~isempty(points_on_slice)
    % 关键：坐标转换
    % 矩阵索引 V(Dim1, Dim2, Dim3) -> 切片 img(Dim2, Dim3)
    % imagesc 绘图时：X轴对应 Dim3，Y轴对应 Dim2
    
    y_coords = points_on_slice(:, 2); % Dim2 (行) -> Y轴
    x_coords = points_on_slice(:, 3); % Dim3 (列) -> X轴
    
    % 绘制实心圆点
    plot(x_coords, y_coords, 'ro', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', ... % 红色实心
        'LineWidth', 1.5);
    
    % (可选) 添加文字标签
    for i = 1:size(points_on_slice, 1)
        text(x_coords(i)+2, y_coords(i), sprintf('Pt%d', i), ...
            'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
    end
else
    warning('没有点位于第 %d 层', target_slice_idx);
end

%% 4. 美化
title(sprintf('Visualization of Slice %d with Marked Points', target_slice_idx));
xlabel('Dimension 3 (Page/Column)');
ylabel('Dimension 2 (Row)');
set(gca, 'TickDir', 'out');
hold off;