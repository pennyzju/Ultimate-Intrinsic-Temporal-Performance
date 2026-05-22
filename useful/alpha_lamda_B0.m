%% 多场强 CV 收敛对比曲线 (1.5T, 3T, 7T)
clear;
B0_list = [1.5, 3, 7];
subdir = '20251022seq';
flag = 'left';
seq_list = [100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200,2500];
%slice_idx = 48;

% 颜色方案 (1.5T: 蓝色, 3T: 绿色, 7T: 红色)
colors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];

% 预分配存储矩阵 [场强数量 x 序列数量]
all_alpha = zeros(numel(B0_list), numel(seq_list));
all_lamda = zeros(numel(B0_list), numel(seq_list));

%% 1. 数据收集阶段
for b = 1:numel(B0_list)
    curr_B0 = B0_list(b);
    % 根据场强动态获取路径
    loss_dir = fullfile(B0_folderpath(curr_B0), subdir, 'loss_plain');
    B1_dir = fullfile(B0_folderpath(curr_B0), subdir, 'B1_map');
    
    fprintf('Processing B0 = %.1f T...\n', curr_B0);
    
    for i = 1:numel(seq_list)
        curr_seq = seq_list(i);
        
        % --- 处理 Alpha 数据 ---
        alpha_file = fullfile(loss_dir, sprintf('alpha_%s_%d.mat', flag, curr_seq));
        if exist(alpha_file, 'file')
            a_data = load(alpha_file);
            if isfield(a_data, 'alpha'), a_val = a_data.alpha; else, a_val = a_data.lamda; end
            temp_a = abs(squeeze(a_val(:, :, :)));
            all_alpha(b, i) = mean(temp_a(temp_a > 0 & isfinite(temp_a)), 'all');
        else
            all_alpha(b, i) = NaN;
        end
        
        % --- 处理 Lamda 数据 ---
        lamda_file = fullfile(B1_dir, sprintf('lamda_%s_%d.mat', flag, curr_seq));
        if exist(lamda_file, 'file')
            l_data = load(lamda_file);
            if isfield(l_data, 'lamda'), l_val = l_data.lamda; else, l_val = l_data.alpha; end
            temp_l = abs(squeeze(l_val(:, :, :)));
            all_lamda(b, i) = mean(temp_l(temp_l > 0 & isfinite(temp_l)), 'all');
        else
            all_lamda(b, i) = NaN;
        end
    end
end

%% 2. 绘图阶段
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 子图 1: Alpha (Loss CV) 收敛对比 ---
ax1 = nexttile; hold on;
for b = 1:numel(B0_list)
    plot(seq_list, all_alpha(b,:), '-o', 'Color', colors(b,:), 'MarkerFaceColor', colors(b,:), ...
        'LineWidth', 2, 'DisplayName', [num2str(B0_list(b)), ' T']);
end
grid on;
title('\alpha (Loss CV) Convergence');
xlabel('Sequence Index (seq)'); ylabel('CV Value');
legend('Location', 'northeast');

% --- 子图 2: Lamda (B1 CV) 收敛对比 ---
ax2 = nexttile; hold on;
for b = 1:numel(B0_list)
    plot(seq_list, all_lamda(b,:), '-s', 'Color', colors(b,:), 'MarkerFaceColor', colors(b,:), ...
        'LineWidth', 2, 'DisplayName', [num2str(B0_list(b)), ' T']);
end
grid on;
title('\lambda (B_1 CV) Convergence');
xlabel('Sequence Index (seq)'); ylabel('CV Value');
legend('Location', 'northeast');

% 统一设置坐标轴刻度字体
set([ax1, ax2], 'FontSize', 11, 'LineWidth', 1);

% 总标题
title(t, sprintf('Convergence Comparison across Field Strengths (Slice %d)', slice_idx), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('Processing Complete.\n');