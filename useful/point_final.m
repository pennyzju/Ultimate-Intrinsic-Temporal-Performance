%% 1. 数据准备
% 请确保路径正确
load('H:\UISNR\20240801_UISNR_output\left_0\default_material_maps.mat')
V = epsilon_r; 
% 为了演示，如果您没有加载数据，这里生成一个随机数据代替
if ~exist('V', 'var')
    V = rand(100, 100, 100) * 80; 
end

target_slice_idx = 48;

% 原始坐标 [Row(层), Col(行/Y), Page(列/X)]
% 对应关系: Center(红), Intermediate(绿), Edge(蓝)
point_list = [
    48 56 80;   % Center Location
    48 70 80;   % Intermediate Location (根据您的代码选的中间点)
    48 83 80    % Edge Location
];

%% 2. 提取切片并旋转
% (1) 提取原始切片
slice_raw = squeeze(V(target_slice_idx, :, :));

% (2) 向左旋转 90 度 (逆时针)
slice_img = rot90(slice_raw); 

%% 3. 坐标变换
points_on_slice = point_list(point_list(:,1) == target_slice_idx, :);

new_x = [];
new_y = [];

if ~isempty(points_on_slice)
    % 获取原始图像的宽度 (用于计算坐标变换)
    old_width = size(slice_raw, 2); 
    
    % 提取原始的绘图坐标
    old_y = points_on_slice(:, 2); % 原始行索引 (Y)
    old_x = points_on_slice(:, 3); % 原始列索引 (X)
    
    % === 坐标变换公式 (逆时针 90 度) ===
    % 新的 X = 原始的 Y
    % 新的 Y = 原始宽度 - 原始 X + 1
    new_x = old_y;
    new_y = old_width - old_x + 1;
end

%% 
% 1. 设置画布：稍微加宽宽度 (由400改为550)，高度不变
% 这样可以让横向排列的图例有呼吸空间
figure('Color', 'w', 'Position', [100, 100, 550, 450]); 

% --- 绘制解剖结构 ---
slice_display = slice_img; 
imagesc(slice_display);
axis image; 
axis off; 

% 保持黑底风格
colormap(gray);     
set(gca, 'Color', 'k'); 

% --- 绘制标记点 ---
hold on;

% 颜色保持：Center(蓝) -> Intermediate(绿) -> Edge(红)
point_colors = {'#1E90FF', '#00FF00', '#FF0000'}; 

% [核心修改]：在标签文字后面加上空格，强制增加间距
point_labels = {'Center        ', 'Intermediate        ', 'Edge'};

for i = 1:length(new_x)
    plot(new_x(i), new_y(i), 'o', ...
        'MarkerSize', 7, ...           
        'MarkerFaceColor', point_colors{i}, ...
        'MarkerEdgeColor', 'w', ...    
        'LineWidth', 1.2);             
end

% --- 图例设置 ---
[leg, icons] = legend(point_labels, ...
    'Location', 'southoutside', ... % 放在图片正下方外侧
    'Orientation', 'horizontal', ...% 水平排列
    'Box', 'off', ...               % 去除边框
    'FontSize', 12, ...             % 字体大小
    'TextColor', 'k');              % 黑色字体

% [微调] 控制图例图标(点)的长度，使其紧凑，把空间留给文字间距
leg.ItemTokenSize = [15, 18]; 

hold off;