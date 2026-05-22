%% ================= Alpha Without Coil Loss (仅保留无损耗) 绘图 =================
clear; clc; 

% 定义基础路径 (如果 B0_folderpath_new 函数不可用，将使用这些路径)
base_paths = {
    'F:\20251201_1p5T_UISNR_output', ...
    'G:\20251201_3T_UISNR_output', ...
    'H:\UISNR\20240801_UISNR_output'
};

B0_list = [1.5, 3, 7];
flags = {'rot', 'left', 'front'}; % 对应 Yaw(Axial), Pitch(Sagittal), Roll(Coronal)
seq_idx = 2500;

% [修改] 仅保留 alpha_wo 的文件名模式
file_pattern = 'alpha_wo_%s_%d.mat'; 

% 定义显示标签
field_titles = {'1.5 T', '3 T', '7 T'};
row_labels = {'Yaw', 'Pitch', 'Roll'}; 

% [关键] 全局色标范围 (根据需要调整，例如 [0, 0.5])
global_clim = [0, 0.2]; 

%% ================= 2. 绘图循环 =================
% [修改] 调整窗口宽度，因为列数变少了 (3列)
figure('Color', 'w', 'Position', [100, 100, 750, 650]); 

% [修改] 建立 3行 x 3列 的布局 (原为 3x6)
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:3 % 遍历行 (视角: Yaw, Pitch, Roll)
    curr_flag = flags{r};
    
    % [修改] 内部循环仅遍历 3 列 (对应 3 个场强)
    for c = 1:3
        p = c; % p 直接对应 B0_list 的索引 (1, 2, 3)
        
        % 1. 构建路径
        try
            % 尝试调用您的自定义函数
            base_path = B0_folderpath_new(B0_list(p)); 
        catch
            % 如果函数不存在，使用上方定义的 base_paths 数组
            if p <= length(base_paths)
                base_path = base_paths{p};
            else
                error('无法确定路径，请检查 B0_folderpath_new 函数或 base_paths');
            end
        end
        
        % [修改] 构建文件名，不再需要判断 if c > 3
        filename = sprintf(file_pattern, curr_flag, seq_idx);
        full_path = fullfile(base_path, '20251022seq', 'loss_plain', filename);
        
        % 2. 读取数据
        img = [];
        if exist(full_path, 'file')
            data = load(full_path);
            
            % 尝试读取变量 (兼容 alpha 或 lamda)
            if isfield(data, 'alpha')
                raw = abs(data.alpha);
            elseif isfield(data, 'lamda')
                raw = abs(data.lamda);
            else
                raw = [];
            end
            
            if ~isempty(raw)
                % 3. 切片提取
                switch curr_flag
                    case 'rot'   % Axial
                        slice = squeeze(raw(4:90, 18:85, 87));
                    case 'left'  % Sagittal
                        slice = squeeze(raw(48, 14:100, 13:104));
                    case 'front' % Coronal
                        slice = squeeze(raw(4:90, 56, 13:104));
                end
                
                % 背景挖空处理
                slice_double = double(slice);
                slice_double(slice_double < 1e-6) = NaN; % 背景设为NaN
                img = rot90(slice_double);
            end
        end
        
        % 4. 绘图
        ax = nexttile;
        if ~isempty(img)
            imagesc(img);
            axis image off;
            
            % 样式设置
            colormap(ax, 'turbo'); 
            clim(global_clim);
            set(ax, 'Color', 'k'); % 设置 NaN 背景为黑色
            
            % [标签] 仅在第一行显示场强标题 (1.5T, 3T, 7T)
            if r == 1
                title(field_titles{p}, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            % [标签] 仅在第一列显示视角名称 (Yaw, Pitch, Roll)
            if c == 1
                ylabel(row_labels{r}, 'Visible', 'on', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', 'k');
            end
        else
            axis off;
        end
    end
end

%% ================= 3. 添加装饰元素 =================

% A. 添加 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\alpha^* (Intrinsic Thermal Noise Variability)';
cb.FontSize = 11;

% B. 添加总标题 (直接加在 Layout 上，简洁明了)
title(t, 'Intrinsic Noise Variability', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', '#0072BD');

% 注意：原代码中的分割线 annotation('line'...) 和左右标题 annotation('textbox'...) 已全部删除