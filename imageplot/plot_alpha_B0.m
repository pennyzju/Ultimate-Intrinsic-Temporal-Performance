%% 多场强 & 多角度 Alpha (Noise Variability) 可视化
clear; clc; close all;

%% ================= 1. 参数设置 =================
% 场强列表 (列)
B0_list = [1.5, 3, 7];

% Flag 列表 (行)
flag_list = {'rot', 'left', 'front'}; 
flag_titles = {'Rotated (Axial)', 'Left (Sagittal)', 'Frontal (Coronal)'};

% 固定参数
seq_idx = 2500;             % 基函数数量
subdir  = '20251022seq';    % 子目录

% 绘图设置
my_colormap = 'jet';        % 颜色表
rotate_img = true;          % 是否旋转图像 (rot90) - 通常 MATLAB 矩阵绘图需要旋转一下才符合人体解剖视角

%% ================= 2. 绘图循环 =================
figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Name', 'Alpha Visualization 3x3');
t = tiledlayout(numel(flag_list), numel(B0_list), 'TileSpacing', 'compact', 'Padding', 'compact');

% 预分配一个 Cell 数组存数据
img_cache = cell(numel(flag_list), numel(B0_list));

% --- 第一遍循环：读取数据 & 动态切片 ---
fprintf('Loading data...\n');
for r = 1:numel(flag_list)     % 遍历行 (Flag)
    curr_flag = flag_list{r};
    
    for c = 1:numel(B0_list)   % 遍历列 (B0)
        curr_B0 = B0_list(c);
        
        % 构建文件路径
        try
            base_path = B0_folderpath_new(curr_B0); 
        catch
            error('请确保 B0_folderpath 函数在路径中，或手动修改 base_path');
        end
        
        file_path = fullfile(base_path, subdir, 'loss_plain', ...
            sprintf('alpha_%s_%d.mat', curr_flag, seq_idx));
        
        % 读取数据
        if exist(file_path, 'file')
            tmp = load(file_path);
            % 优先获取 alpha，否则 lamda
            if isfield(tmp, 'alpha')
                data_vol = abs(tmp.alpha);
            elseif isfield(tmp, 'lamda')
                data_vol = abs(tmp.lamda);
            else
                data_vol = [];
            end
            
            % === 【关键修改】根据 Flag 动态截取不同层面和范围 ===
            if ~isempty(data_vol)
                sl = []; % 初始化
                
                switch curr_flag
                    case 'left'
                        % 需求: (48, 14:100, 13:104) -> 侧面 Sagittal
                        % 维度1固定，取维度2和3的特定范围
                        try
                            sl = squeeze(data_vol(48, 14:100, 13:104));
                        catch
                            warning('Left 索引越界，使用默认中间层');
                            sl = squeeze(data_vol(48, :, :));
                        end
                        
                    case 'front'
                        % 需求: (4:90, 56, 13:104) -> 正面 Coronal
                        % 维度2固定，取维度1和3的特定范围
                        try
                            sl = squeeze(data_vol(4:90, 56, 13:104));
                        catch
                            warning('Front 索引越界，使用默认中间层');
                            sl = squeeze(data_vol(:, 56, :));
                        end
                        
                    case {'rot', 'rot_new'} % 兼容可能的命名
                        % 需求: (4:90, 13:94, 87) -> 横断 Axial (假设 rot 是 Z轴方向)
                        % 维度3固定，取维度1和2的特定范围
                        try
                            sl = squeeze(data_vol(4:90, 13:94, 87));
                        catch
                            warning('Rot 索引越界，使用默认中间层');
                            sl = squeeze(data_vol(:, :, 87));
                        end
                        
                    otherwise
                        warning('未知的 flag: %s', curr_flag);
                end
                
                % 旋转处理 (根据需要调整，squeeze 后的矩阵方向通常需要转一下)
                if rotate_img
                    sl = rot90(sl); 
                end
                
                img_cache{r, c} = sl;
            end
        else
            fprintf('File not found: %s\n', file_path);
        end
    end
end

% --- 第二遍循环：绘图 (按行统一 Clim) ---
fprintf('Plotting...\n');
for r = 1:numel(flag_list)
    
    % 计算当前行 (同一位置，不同场强) 的 Clim
    row_data = [img_cache{r, :}]; 
    if isempty(row_data)
        clim_range = [0 1];
    else
        % 使用 99% 分位数自动去除噪点极值，使对比度更好
        clim_range = [0, prctile(row_data(:), 99)]; 
    end
    
    for c = 1:numel(B0_list)
        nexttile(t);
        img = img_cache{r, c};
        
        if ~isempty(img)
            imagesc(img);
            axis image; axis off;
            colormap(my_colormap);
            clim(clim_range); % 应用统一色标
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
        end
        
        % 第一行显示场强
        if r == 1
            title(sprintf('%.1f T', B0_list(c)), 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        % 第一列显示 Flag 名称
        if c == 1
            % 使用 text 在左侧显示，比 ylabel 更灵活
            text(-0.1, 0.5, flag_titles{r}, 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'Rotation', 90, ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
        
        % 最后一列显示 Colorbar
        if c == numel(B0_list)
            cb = colorbar;
            cb.Label.String = '\alpha^*';
        end
    end
end

% 总标题
title(t, sprintf('Alpha Map Comparison (Seq=%d)', seq_idx), 'FontSize', 14, 'FontWeight', 'bold');