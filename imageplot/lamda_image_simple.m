clear; clc;

%% ================= 1. 全局配置 (修改这里即可) =================
% [核心参数] 手动设置统一色标范围
% 建议设为 [0, 0.12] 或 [0, 0.15]。
% 设得越低，1.5T/3T 越亮，但 7T 的红区会更饱和。
%global_clim = [-10, -1]; 
global_clim = [0 0.15]; %线性范围，修改64行

% 路径配置 (顺序: 1.5T, 3T, 7T)
base_paths = {
    'F:\20251201_1p5T_UISNR_output', ...
    'G:\20251201_3T_UISNR_output', ...
    'H:\UISNR\20240801_UISNR_output'
};

% 维度配置 (行顺序: Yaw/Rot, Pitch/Left, Roll/Front)
flags = {'rot', 'left', 'front'}; 

% 文件名模板 (对应左侧3列和右侧3列)
file_patterns = {'lamda_%s_2500.mat'};

%% ================= 2. 自动化绘图循环 (SCI 优化版) =================
figure('Color', 'w', 'Position', [100, 100, 600, 600]); % 画布稍微大一点
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');



% 定义标题和标签以便循环使用
field_titles = {'1.5 T', '3 T', '7 T'};
row_labels = {'Yaw', 'Pitch','Roll' }; % 对应 rot(Z), left(X), front(Y)

for r = 1:3 % 遍历行 (Flag: rot, left, front)
    curr_flag = flags{r};
    
    curr_pattern = file_patterns{1};
        
    for p = 1:3 % 遍历场强: 1.5T, 3T, 7T
        % 1. 构建文件路径
        full_path = fullfile(base_paths{p}, '20251022seq', 'B1_map', ...
            sprintf(curr_pattern, curr_flag));
        
        % 2. 读取数据 (带错误保护)
        img = [];
        if exist(full_path, 'file')
            data = load(full_path);
            if isfield(data, 'lamda')
                % 3. 根据 flag 自动选择切片逻辑
                raw = data.lamda;
                switch curr_flag
                    case 'rot'   % Axial
                        slice = squeeze(raw(4:90, 18:85, 87));
                    case 'left'  % Sagittal
                        slice = squeeze(raw(48, 14:100, 13:104));
                    case 'front' % Coronal
                        slice = squeeze(raw(4:90, 56, 13:104));
                end
                
                % [关键修改 A] 背景挖空：将接近0的值设为 NaN
                % 这样它们就会透明，显示出底部的颜色
                slice_double = double(slice);
                slice_double(slice_double < 1e-6) = NaN; 
                
                img = rot90(slice_double); 
            end
        end
        
        % 4. 绘图 (自动分配位置)
        ax = nexttile; 
        if ~isempty(img)
            imagesc(img);
            axis image off;
            
            % [关键修改 B] 设置色阶和范围
            colormap(ax, 'turbo'); % 使用 turbo 增强对比
            clim(global_clim);     % 强制应用统一范围
            
            % [关键修改 C] 设置 Axes 背景为黑色
            % 这样 NaN 区域就会显示为黑色
            set(ax, 'Color', 'k'); 
            
            % [美化] 仅在第一行添加场强标题
            if r == 1
                title(field_titles{p}, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            % [美化] 仅在第一列添加左侧视角标签
            if p == 1
                % 使用 text 或 ylabel 技巧
                ylabel(row_labels{r}, 'Visible', 'on', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
                % 如果 ylabel 被 axis off 隐藏了，可以用 text 手动加：
                % text(-10, size(img,1)/2, row_labels{r}, ...
                %    'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Rotation', 90);
            end
        else
            % 如果文件缺失，画个空图占位
            axis off; 
        end
    end
end

%% ================= 3. 一步添加共享 Colorbar =================
cb = colorbar; 
cb.Layout.Tile = 'east'; % 将 Colorbar 放置在布局的最右侧
cb.Label.String = '\lambda^* (Intrinsic Sensitivity Variability)'; % 添加物理量标签
cb.Label.FontSize = 12;
cb.FontSize = 11;