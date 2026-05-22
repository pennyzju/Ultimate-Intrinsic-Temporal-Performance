%% 多场强 & 多角度 Alpha (Noise Variability) 可视化
clear; clc; close all;

%% ================= 1. 参数设置 =================
% 场强列表 (列)
B0_list = [1.5, 3, 7];

% Flag 列表 (行) - 对应文件名中的 flag 部分
flag_list = {'rot', 'left', 'front'}; 
% 对应的显示标题 (可选)
flag_titles = {'Rotated', 'Left', 'Frontal'};

% 固定参数
seq_idx = 2500;             % 基函数数量
subdir  = '20251022seq';    % 子目录
target_slice = 48;          % 切片索引 (根据您的模型调整，例如 48)
slice_dim = 1;              % 切片维度: 1=Sagittal(侧面), 2=Coronal(冠状), 3=Axial(横断)
                            % 注意：根据您之前的代码，数据维度可能是 [Row, Col, Page]

% 绘图设置
my_colormap = 'jet';        % 颜色表
rotate_img = true;          % 是否旋转图像 (rot90)

%% ================= 2. 绘图循环 =================
figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Name', 'Alpha Visualization 3x3');
t = tiledlayout(numel(flag_list), numel(B0_list), 'TileSpacing', 'compact', 'Padding', 'compact');

% 预分配一个 Cell 数组存数据，以便后续统一 Clim
img_cache = cell(numel(flag_list), numel(B0_list));

% --- 第一遍循环：读取数据 ---
fprintf('Loading data...\n');
for r = 1:numel(flag_list)     % 遍历行 (Flag)
    curr_flag = flag_list{r};
    
    for c = 1:numel(B0_list)   % 遍历列 (B0)
        curr_B0 = B0_list(c);
        
        % 构建文件路径
        % 假设 B0_folderpath 是您的自定义函数
        try
            base_path = B0_folderpath_new(curr_B0); 
        catch
            error('请确保 B0_folderpath 函数在路径中，或手动修改 base_path');
        end
        
        file_path = fullfile(base_path, subdir, 'loss_plain', ...
            sprintf('alpha_wo_%s_%d.mat', curr_flag, seq_idx));
        
        % 读取数据
        if exist(file_path, 'file')
            tmp = load(file_path);
            % 尝试获取变量，优先 alpha，其次尝试 lamda (防错)
            if isfield(tmp, 'alpha')
                data_vol = abs(tmp.alpha);
            elseif isfield(tmp, 'lamda')
                data_vol = abs(tmp.lamda);
                fprintf('Warning: Variable "alpha" not found in %s, used "lamda" instead.\n', file_path);
            else
                warning('No valid variable found in %s', file_path);
                data_vol = [];
            end
            
            % 提取切片
            if ~isempty(data_vol)
                switch slice_dim
                    case 1
                        sl = squeeze(data_vol(target_slice, :, :));
                    case 2
                        sl = squeeze(data_vol(:, target_slice, :));
                    case 3
                        sl = squeeze(data_vol(:, :, target_slice));
                end
                
                if rotate_img
                    sl = rot90(sl); % 根据需要旋转
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
    
    % 计算当前行 (同一位置，不同场强) 的最大值，用于统一色标
    row_data = [img_cache{r, :}]; % 拼接这一行所有数据
    if isempty(row_data)
        clim_range = [0 1];
    else
        % 取 98% 分位数作为上限，避免极值（热点）掩盖整体细节
        % 或者直接用 max(row_data(:))
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
        
        % 仅在第一行显示场强标题
        if r == 1
            title(sprintf('%.1f T', B0_list(c)), 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        % 仅在第一列显示 Flag 标题 (Y轴标签效果)
        if c == 1
            ylabel(flag_titles{r}, 'FontSize', 12, 'FontWeight', 'bold', 'Visible', 'on');
            % 去除刻度只留 Label 可能会有点麻烦，直接用 text 写在左边也可以
            % 这里简单处理：不做额外 text，直接利用 ylabel
        end
        
        % 仅在最后一列显示 Colorbar
        if c == numel(B0_list)
            cb = colorbar;
            cb.Label.String = '\alpha^* (Noise Var.)';
        end
    end
end

% 添加总标题
title(t, sprintf('Intrinsic Noise Variability (\\alpha^*) Comparison (Seq=%d)', seq_idx), ...
    'FontSize', 14, 'FontWeight', 'bold');