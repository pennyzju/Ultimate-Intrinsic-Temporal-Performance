
%% ================= Alpha 包含损耗 vs 不含损耗 对比绘图 =================
clear; clc; 

base_paths = {
    'F:\20251201_1p5T_UISNR_output', ...
    'G:\20251201_3T_UISNR_output', ...
    'H:\UISNR\20240801_UISNR_output'
};
B0_list = [1.5, 3, 7];
flags = {'rot', 'left', 'front'}; % 对应 Axial, Sagittal, Coronal
seq_idx = 2500;

% 定义文件名模式
% 假设: pattern{1} 是无损耗 (alpha_wo), pattern{2} 是有损耗 (alpha)
file_patterns = {'alpha_wo_%s_%d.mat', 'alpha_%s_%d.mat'}; 

% 定义显示标签
field_titles = {'1.5 T', '3 T', '7 T'};
% 注意：根据之前的切片逻辑，rot通常是Axial，left是Sagittal，front是Coronal
row_labels = {'Yaw', 'Pitch', 'Roll'}; 

% [关键] 全局色标范围
% Alpha 通常数值较小，建议设为 [0, 0.05] 或 [0, 0.1] 以便看清差异
% 如果有损耗(右边)的值很大，可能需要适当调高上限
%global_clim = [0, 0.1]; 
global_clim = [0,0.5]; 
%% ================= 2. 绘图循环 =================
figure('Color', 'w', 'Position', [50, 50, 1400, 700]);
% 建立 3行 x 6列 的布局
t = tiledlayout(3, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:3 % 遍历行 (视角: Axial, Sagittal, Coronal)
    curr_flag = flags{r};
    
    % 内部循环遍历 6 列
    for c = 1:6
        % --- 逻辑判断：当前列属于哪种类型，哪个场强 ---
        if c <= 3
            % 左边 3 列：Without Coil Loss
            type_idx = 1; 
            p = c; % p = 1, 2, 3
        else
            % 右边 3 列：With Coil Loss
            type_idx = 2;
            p = c - 3; % p = 1, 2, 3
        end
        
        curr_pattern = file_patterns{type_idx};
        
        % 1. 构建路径 (请确保 base_paths 函数可用)
        try
            base_path = B0_folderpath_new(B0_list(p)); 
        catch
            error('请定义 B0_folderpath 或手动设置路径');
        end
        
        full_path = fullfile(base_path, '20251022seq', 'loss_plain', ...
            sprintf(curr_pattern, curr_flag, seq_idx));
        
        % 2. 读取数据
        img = [];
        if exist(full_path, 'file')
            data = load(full_path);
            % 尝试读取 alpha (兼容不同变量名)
            if isfield(data, 'alpha')
                raw = abs(data.alpha);
            elseif isfield(data, 'lamda')
                raw = abs(data.lamda);
            else
                raw = [];
            end
            
            if ~isempty(raw)
                % 3. 切片提取 (保持原逻辑)
                switch curr_flag
                    case 'rot'   % Axial
                        slice = squeeze(raw(4:90, 18:85, 87));
                    case 'left'  % Sagittal
                        slice = squeeze(raw(48, 14:100, 13:104));
                    case 'front' % Coronal
                        slice = squeeze(raw(4:90, 56, 13:104));
                end
                
                % [关键] 背景挖空处理 (黑底效果)
                slice_double = double(slice);
                slice_double(slice_double < 1e-6) = NaN; % 背景设为NaN
                img = rot90(slice_double);%这里设置是否为log
            end
        end
        
        % 4. 绘图
        ax = nexttile;
        if ~isempty(img)
            imagesc(img);
            axis image off;
            
            % 样式设置
            colormap(ax, 'turbo'); % 或者 'turbo'
            clim(global_clim);
            set(ax, 'Color', 'k'); % 设置 NaN 背景为黑色
            
            % [标签] 仅在第一行显示场强
            if r == 1
                title(field_titles{p}, 'FontSize', 12, 'FontWeight', 'bold');
            end
            
            % [标签] 仅在第一列显示视角名称
            if c == 1
                ylabel(row_labels{r}, 'Visible', 'on', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', 'k');
            end
        else
            axis off;
        end
    end
end

%% ================= 3. 添加装饰元素 (分割线与大标题) =================

% A. 添加共享 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\alpha^* (Intrinsic Thermal Noise Variability)';
cb.FontSize = 11;

% B. 绘制中间分割线 (模拟图中虚线)
% 获取整个 layout 的像素位置，计算中间位置
% 这里使用 annotation 绘制一条竖线
% [x_start, y_start, x_width, y_height] (归一化坐标 0-1)
annotation('line', [0.495 0.495], [0.05 0.92], ...
    'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);

% C. 添加顶部大标题 (Without / With Coil Loss)
% 左侧标题
annotation('textbox', [0.15, 0.94, 0.3, 0.05], ...
    'String', 'Without Coil Loss', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', '#0072BD'); % 蓝色

% 右侧标题
annotation('textbox', [0.58, 0.94, 0.3, 0.05], ...
    'String', 'With Coil Loss', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', '#D95319'); % 橙色/红色