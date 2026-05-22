%% 1. 数据准备
load('H:\UISNR\20240801_UISNR_output\left_0\default_material_maps.mat')
V = epsilon_r;

target_slice_idx = 48;

% 原始坐标 [Row(层), Col(行/Y), Page(列/X)]
point_list = [
    48 56 80;   % Center Location
    48 67 80;   % Intermediate Location (根据您的代码选的中间点)
    48 75 80    % Edge Location
];

%% 2. 提取切片并旋转
% (1) 提取原始切片
slice_raw = squeeze(V(target_slice_idx, :, :));

% (2) 向左旋转 90 度 (逆时针)
% rot90(A) 默认就是逆时针 90 度
% 如果想向右转(顺时针)，用 rot90(A, -1)
slice_img = rot90(slice_raw); 

%% 3. 坐标变换 (核心步骤)
% 筛选出当前层的点
points_on_slice = point_list(point_list(:,1) == target_slice_idx, :);

if ~isempty(points_on_slice)
    % 获取原始图像的宽度 (用于计算坐标变换)
    % 注意：要在旋转前获取尺寸，或者旋转后取高度
    old_width = size(slice_raw, 2); 
    
    % 提取原始的绘图坐标 (注意 imagesc 的 x,y 对应关系)
    old_y = points_on_slice(:, 2); % 原始行索引 (Y)
    old_x = points_on_slice(:, 3); % 原始列索引 (X)
    
    % === 坐标变换公式 (逆时针 90 度) ===
    % 新的 X = 原始的 Y
    % 新的 Y = 原始宽度 - 原始 X + 1
    new_x = old_y;
    new_y = old_width - old_x + 1;
end

%% 4. 可视化
figure('Color', 'w');

% 绘制旋转后的图像
imagesc(slice_img);
colormap('gray');
axis image; % 这一步很重要，保证旋转后比例不拉伸
colorbar;

hold on;

if ~isempty(points_on_slice)
    % 使用变换后的新坐标绘图
    plot(new_x, new_y, 'ro', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', ...
        'LineWidth', 1.5);
    
    % 标记文字
    for i = 1:length(new_x)
        text(new_x(i)+2, new_y(i), sprintf('Pt%d', i), ...
            'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

%title(sprintf('Slice %d (Rotated 90^o Left)', target_slice_idx));
% 注意：旋转后，坐标轴含义变了
xlabel('Original Row Dimension');
ylabel('Original Column Dimension (Inverted)');
set(gca, 'TickDir', 'out');
hold off;