clear;
foldername = 'left_0';

% --- 1. 加载数据 ---
basePath = 'H:\UISNR\20240801_UISNR_output';
maskFilename = fullfile(basePath, foldername, [foldername '_material_maps.mat']);
maskdipFilename = fullfile(basePath, 'public/default_material_maps.mat');

% 加载
data1 = load(maskFilename, 'rotated_epsilon_r');
rotated_mask = data1.rotated_epsilon_r;
data2 = load(maskdipFilename, 'mask_dip2');
mask_dip2 = data2.mask_dip2;

% --- 2. 提取中间层 ---
[d1, d2, d3] = size(rotated_mask);
slice_dim = 1; % 侧视图
mid_idx = round(d1 / 2);

% 提取切片
slice_rot = squeeze(rotated_mask(mid_idx, :, :));
slice_dip = squeeze(mask_dip2(mid_idx, :, :));

% 旋转90度以符合视觉
slice_rot = rot90(slice_rot);
slice_dip = rot90(slice_dip);

% --- 3. 准备顶层：Standard Mask (深红色厚轮廓) ---
mask_bin = slice_dip > 0;
perim = bwperim(mask_bin);

% 制造"厚度"
thickness_radius = 2; 
se = strel('disk', thickness_radius);
thick_band = imdilate(perim, se);

% 设置暗红色覆盖层 (R通道=0.6)
overlay_color = zeros([size(slice_rot), 3]); 
overlay_color(:, :, 1) = 0.6;  

% --- 4. 绘图与叠加 ---
figure('Color', 'k', 'Position', [100, 100, 600, 600]);

% 层 1: 绘制内部连续数值矩阵
imagesc(slice_rot);

% 修改色图，让背景值(1)完美融合黑底
cmap = jet(256); 
cmap(1, :) = [0 0 0]; 
colormap(gca, cmap);

% 锁定颜色范围
bg_val = min(slice_rot(:)); 
clim([bg_val, max(slice_rot(:))]); 

% % 添加 Colorbar 
% cb = colorbar;
% cb.Color = 'w'; 

axis image off;
hold on;

% 层 2: 绘制"暗红厚壳"
% h_overlay = image(overlay_color);
% alpha_map = double(thick_band) * 0.9; 
% set(h_overlay, 'AlphaData', alpha_map);

% --- 5. 标记指定的 3D 数据点 ---
point_list = [
    48 56 80;   % Center Location
    48 70 80;   % Intermediate Location (根据您的代码选的中间点)
    48 83 80    % Edge Location
];

% 坐标映射
pts_y = point_list(:, 2);
pts_z = point_list(:, 3);
plot_X = pts_y;
plot_Y = d3 - pts_z + 1;

% 定义形状与标签
marker_shapes = {'^', 's', 'o'};
marker_labels = {'Center', 'Intermediate', 'Edge'};

% 清除 hold 状态下的旧句柄，确保重新开始
hold on; 

% 循环绘制并直接存储句柄
h_pts = []; % 初始化为空数组
for i = 1:length(marker_shapes)
    m_size = 150;
    if strcmp(marker_shapes{i}, 'p'), m_size = 140; end
    
    % 绘制点并保存句柄到数组
    h_pts(i) = scatter(plot_X(i), plot_Y(i), m_size, ...
        'Marker', marker_shapes{i}, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'w', ...
        'DisplayName', marker_labels{i}); % 直接在 scatter 中绑定标签
end

% --- 6. 添加图例 ---
% 使用绑定好的 'DisplayName' 自动生成图例，防止句柄失效
lgd = legend(h_pts); 

% 设置图例样式
set(lgd, ...
    'Color', 'none', ...         % 背景透明
    'TextColor', 'w', ...        % 文字白色
    'EdgeColor', 'w', ...        % 框线白色
    'FontSize', 9, ...
    'Location', 'northeast');
% % --- 6. 装饰与保存 ---
% title({['Comparison Slice: ' num2str(mid_idx)], ...
%        'Rotated Anatomy vs. Boundary with 3 Data Points'}, ...
%        'Color', 'w', 'FontSize', 12);

outfile = sprintf('%s_overlay_thick_colored_points.png', foldername);
fprintf('Saved overlay visualization to: %s\n', outfile);
saveas(gcf, outfile);







