clear;
foldername = 'rot_-5';

% ==========================================
% 【用户配置区】
% 在这里选择你想看的视角: 'left', 'front', 或 'rot'
view_mode = 'rot';  
% ==========================================

% --- 1. 加载数据 ---
basePath = 'H:\UISNR\20240801_UISNR_output';
maskFilename = fullfile(basePath, foldername, [foldername '_material_maps.mat']);
maskdipFilename = fullfile(basePath, 'public/default_material_maps.mat');

% 加载
data1 = load(maskFilename, 'rotated_epsilon_r');
rotated_mask = data1.rotated_epsilon_r;
data2 = load(maskdipFilename, 'mask_dip2');
mask_dip2 = data2.mask_dip2;

% --- 2. 获取矩阵维度与全局最大值 ---
[d1, d2, d3] = size(rotated_mask);
mid_x = round(d1 / 2);
mid_y = round(d2 / 2);
mid_z = round(d3 / 2);

% 获取全局最大值，确保无论切哪个面，Colorbar 颜色标准一致
global_max = max(rotated_mask(:));

% --- 3. 根据选择的视角提取切片 ---
switch lower(view_mode)
    case 'left'
        % 侧视图 (固定 X 轴)
        slice_rot = squeeze(rotated_mask(mid_x, :, :));
        slice_dip = squeeze(mask_dip2(mid_x, :, :));
        current_idx = mid_x;
        view_title = 'Left View (Sagittal)';
        
    case 'front'
        % 正视图 (固定 Y 轴)
        slice_rot = squeeze(rotated_mask(:, mid_y, :));
        slice_dip = squeeze(mask_dip2(:, mid_y, :));
        current_idx = mid_y;
        view_title = 'Front View (Coronal)';
        
    case 'rot'
        % 俯视图/轴状图 (固定 Z 轴)
        slice_rot = squeeze(rotated_mask(:, :, mid_z));
        slice_dip = squeeze(mask_dip2(:, :, mid_z));
        current_idx = mid_z;
        view_title = 'Rot View (Axial)';
        
    otherwise
        error('未知的视角，请将 view_mode 设置为 ''left'', ''front'', 或 ''rot''');
end

% 旋转90度以符合视觉习惯
slice_rot = rot90(slice_rot);
slice_dip = rot90(slice_dip);

% 【核心修复】：将空气及旋转填充产生的 0 统一处理，解决蓝底问题
slice_rot(slice_rot <= 1.01) = 0;

% --- 4. 准备顶层：Standard Mask (深红色厚轮廓) ---
mask_bin = slice_dip > 0;
perim = bwperim(mask_bin);

% 制造"厚度"
thickness_radius = 2; 
se = strel('disk', thickness_radius);
thick_band = imdilate(perim, se);

% 设置暗红色覆盖层 (R通道=0.6)
overlay_color = zeros([size(slice_rot), 3]); 
overlay_color(:, :, 1) = 0.6;  

% --- 5. 绘图与叠加 ---
figure('Color', 'k', 'Position', [100, 100, 600, 600]);

% 绘制内部连续数值矩阵
imagesc(slice_rot);

% 修改色图，将第一阶颜色强制设为纯黑
cmap = jet(256); 
cmap(1, :) = [0 0 0]; 
colormap(gca, cmap);

% 锁定 0 为下限，全局最大值为上限
clim([0, global_max]); 

% % 添加 Colorbar 
% cb = colorbar;
% cb.Color = 'w'; 
% cb.Label.String = 'Epsilon_R';
% cb.Label.Color = 'w';

axis image off;
hold on;

% 绘制"暗红厚壳"
h_overlay = image(overlay_color);
alpha_map = double(thick_band) * 0.9; 
set(h_overlay, 'AlphaData', alpha_map);

% % 装饰标题
% title({sprintf('%s - Slice: %d', view_title, current_idx), ...
%        'Rotated Anatomy vs. Standard Boundary'}, ...
%        'Color', 'w', 'FontSize', 14);

% --- 6. 动态保存 ---
% 根据选择的视角命名保存的文件，避免相互覆盖
outfile = sprintf('%s_overlay_%s.png', foldername, view_mode);
fprintf('Saved visualization to: %s\n', outfile);
saveas(gcf, outfile);