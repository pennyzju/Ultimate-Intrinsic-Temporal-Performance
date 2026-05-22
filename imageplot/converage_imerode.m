%% 多场强 CV 收敛对比曲线 - 指定Slice(48)层面缩小5像素处理
clear;
B0_list = [1.5, 3, 7];
subdir = '20251022seq';
flag = 'left';
seq_list = [100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200, 2500];

% ==================== 关键修改点 ====================
slice_idx = 48;        % 指定有数据的层数
shrink_pixels = 5;     % 向内收缩的像素数
% ====================================================

% 颜色方案
colors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];

% 预分配
all_alpha = zeros(numel(B0_list), numel(seq_list));
all_lamda = zeros(numel(B0_list), numel(seq_list));

%% 1. 数据收集阶段
for b = 1:numel(B0_list)
    curr_B0 = B0_list(b);
    % 请确保 B0_folderpath 函数在您的路径中
    loss_dir = fullfile(B0_folderpath(curr_B0), subdir, 'loss_plain');
    B1_dir = fullfile(B0_folderpath(curr_B0), subdir, 'B1_map');
    
    fprintf('Processing B0 = %.1f T (Slice %d)...\n', curr_B0, slice_idx);
    
    % 初始化当前场强的 Mask
    alpha_mask_2d = [];
    lamda_mask_2d = [];
    
    for i = 1:numel(seq_list)
        curr_seq = seq_list(i);
        
        % ----------------- 处理 Alpha 数据 -----------------
        alpha_file = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, curr_seq));
        if exist(alpha_file, 'file')
            a_data = load(alpha_file);
            if isfield(a_data, 'alpha'), a_val = a_data.alpha; else, a_val = a_data.lamda; end
            
            % 【关键步骤 1】提取指定层面的 2D 数据
            % 假设数据维度是 [X, Y, Z] 或 [Slice, H, W]，提取第48行
            temp_3d = abs(a_val);
            temp_2d = squeeze(temp_3d(slice_idx, :, :)); 
            
            % 【关键步骤 2】生成并腐蚀 2D Mask (仅首次执行)
            if isempty(alpha_mask_2d)
                % 生成初始 Mask：大于0且非无穷
                raw_mask = (temp_2d > 0) & isfinite(temp_2d);
                
                % 定义 2D 结构元素 (3x3 正方形)
                % 每次腐蚀会剥离最外层 1 个像素
                se = strel('square', 3); 
                
                % 循环腐蚀 5 次
                current_mask = raw_mask;
                for k = 1:shrink_pixels
                    current_mask = imerode(current_mask, se);
                end
                alpha_mask_2d = current_mask;
                
                % (可选) 调试显示一下 Mask 也就是看看切得对不对
                % figure; imshow(alpha_mask_2d); title(['Alpha Mask B0=', num2str(curr_B0)]);
            end
            
            % 【关键步骤 3】利用 2D Mask 提取数值并求均值
            valid_vals = temp_2d(alpha_mask_2d);
            
            if ~isempty(valid_vals)
                all_alpha(b, i) = mean(valid_vals, 'all');
            else
                all_alpha(b, i) = NaN;
                fprintf('Warning: Alpha Mask is empty after shrinking at seq %d\n', curr_seq);
            end
        else
            all_alpha(b, i) = NaN;
        end
        
        % ----------------- 处理 Lamda 数据 -----------------
        lamda_file = fullfile(B1_dir, sprintf('lamda_%s_%d.mat', flag, curr_seq));
        if exist(lamda_file, 'file')
            l_data = load(lamda_file);
            if isfield(l_data, 'lamda'), l_val = l_data.lamda; else, l_val = l_data.alpha; end
            
            temp_3d = abs(l_val);
            temp_2d = squeeze(temp_3d(slice_idx, :, :)); % 提取 2D 切片
            
            if isempty(lamda_mask_2d)
                raw_mask = (temp_2d > 0) & isfinite(temp_2d);
                se = strel('square', 3);
                
                current_mask = raw_mask;
                for k = 1:shrink_pixels
                    current_mask = imerode(current_mask, se);
                end
                lamda_mask_2d = current_mask;
            end
            
            valid_vals = temp_2d(lamda_mask_2d);
            
            if ~isempty(valid_vals)
                all_lamda(b, i) = mean(valid_vals, 'all');
            else
                all_lamda(b, i) = NaN;
            end
        else
            all_lamda(b, i) = NaN;
        end
    end
end

%% 2. 绘图阶段
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: Alpha
ax1 = nexttile; hold on;
for b = 1:numel(B0_list)
    plot(seq_list, all_alpha(b,:), '-o', 'Color', colors(b,:), 'MarkerFaceColor', colors(b,:), ...
        'LineWidth', 2, 'DisplayName', [num2str(B0_list(b)), ' T']);
end
grid on;
title(sprintf('\\alpha Convergence (Slice %d, Shrink %d px)', slice_idx, shrink_pixels));
xlabel('Sequence Index'); ylabel('CV Value');
legend('Location', 'northeast');

% 子图 2: Lamda
ax2 = nexttile; hold on;
for b = 1:numel(B0_list)
    plot(seq_list, all_lamda(b,:), '-s', 'Color', colors(b,:), 'MarkerFaceColor', colors(b,:), ...
        'LineWidth', 2, 'DisplayName', [num2str(B0_list(b)), ' T']);
end
grid on;
title(sprintf('\\lambda Convergence (Slice %d, Shrink %d px)', slice_idx, shrink_pixels));
xlabel('Sequence Index'); ylabel('CV Value');
legend('Location', 'northeast');

set([ax1, ax2], 'FontSize', 11, 'LineWidth', 1);
title(t, 'Convergence Comparison on Specific Slice', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Processing Complete.\n');