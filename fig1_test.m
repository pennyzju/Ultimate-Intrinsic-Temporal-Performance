%% 多场强 CV 收敛对比（方案 B：固定点，点级收敛）- SCI 优化版
clear; clc; close all;

%% ================= 基本参数 =================
% 请确保 B0_folderpath 函数在您的路径中可用
B0_list = [1.5, 3, 7];
subdir  = '20251022seq';
flag    = 'left';

% 序列列表
seq_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200, 2500];

%% ================= 方案 B：人工指定代表性点 =================
% 每一行是一个体素坐标 [x y z]
% 对应颜色逻辑: Center(蓝) -> Intermediate(绿) -> Edge(红)
point_list = [
    48 56 74;   % Point 1: Center
    48 70 74;   % Point 2: Intermediate
    48 84 74;   % Point 3: Edge (Left)
];

point_names = {'Center', 'Intermediate', 'Edge'};
% 定义 SCI 绘图标准色 (Hex 格式)
% Center: 亮蓝 (#1E90FF), Intermediate: 亮绿 (#00FF00), Edge: 纯红 (#FF0000)
my_colors = {'#1E90FF', '#00FF00', '#FF0000'}; 

num_points = size(point_list,1);

%% ================= Step 1 & 2：数据收集 (保持原逻辑) =================
% --- 1.1 获取线性索引 ---
test_B0  = B0_list(1);
test_seq = seq_list(1);
try
    test_dir = fullfile(B0_folderpath(test_B0), subdir, 'loss_plain');
catch
    error('请确保 B0_folderpath 函数在 MATLAB 路径中。');
end
test_file = fullfile(test_dir, sprintf('alpha_wo_%s_%d.mat', flag, test_seq));
assert(exist(test_file,'file')==2, 'Reference file not found: %s', test_file);

tmp = load(test_file);
if isfield(tmp,'alpha'), ref = tmp.alpha; else, ref = tmp.lamda; end
data_size = size(ref);

sample_idx = sub2ind(data_size, point_list(:,1), point_list(:,2), point_list(:,3));

% --- 1.2 预分配 ---
all_alpha = nan(numel(B0_list), num_points, numel(seq_list));
all_lamda = nan(numel(B0_list), num_points, numel(seq_list));

% --- 1.3 循环读取 ---
for b = 1:numel(B0_list)
    curr_B0 = B0_list(b);
    loss_dir = fullfile(B0_folderpath(curr_B0), subdir, 'loss_plain');
    B1_dir   = fullfile(B0_folderpath(curr_B0), subdir, 'B1_map');
    
    fprintf('Processing B0 = %.1f T...\n', curr_B0);
    
    for i = 1:numel(seq_list)
        curr_seq = seq_list(i);
        
        % Alpha 读取
        alpha_file = fullfile(loss_dir, sprintf('alpha_%s_%d.mat', flag, curr_seq));
        if exist(alpha_file,'file')
            a_data = load(alpha_file);
            if isfield(a_data,'alpha'), temp_a = abs(a_data.alpha); else, temp_a = abs(a_data.lamda); end
            vals = temp_a(sample_idx);
            vals(~isfinite(vals) | vals <= 0) = NaN;
            all_alpha(b,:,i) = vals;
        end
        
        % Lamda 读取
        lamda_file = fullfile(B1_dir, sprintf('lamda_%s_%d.mat', flag, curr_seq));
        if exist(lamda_file,'file')
            l_data = load(lamda_file);
            if isfield(l_data,'lamda'), temp_l = abs(l_data.lamda); else, temp_l = abs(l_data.alpha); end
            vals = temp_l(sample_idx);
            vals(~isfinite(vals) | vals <= 0) = NaN;
            all_lamda(b,:,i) = vals;
        end
    end
end
fprintf('Data processing complete.\n');

%% ================= Step 3：可视化绘图 (SCI 风格优化) =================

% 创建一个大图 (2行 x 3列)
figure('Color','w','Position',[100 100 1400 700], 'Name', 'Metric Convergence Comparison');
t = tiledlayout(2, numel(B0_list), 'TileSpacing','compact', 'Padding','compact');

% --- 3.1 计算全局统一 Y 轴范围 (便于横向对比) ---
% Lamda 范围
valid_l = all_lamda(~isnan(all_lamda));
ylim_lamda = [0, max(valid_l(:)) * 1.1]; 

% Alpha 范围 (注意：Edge 在 7T 可能很大，如果想看细节，可以用对数坐标或截断)
valid_a = all_alpha(~isnan(all_alpha));
ylim_alpha = [0, max(valid_a(:)) * 1.1];
% 如果 Alpha 在 7T Edge 处极其巨大，建议这里手动 clamp 一下，例如:
% ylim_alpha = [0, 1.2]; 

%% --- 3.2 绘制第一行：Lamda (Sensitivity Variability) ---
for b = 1:numel(B0_list)
    ax = nexttile(t); hold on;
    
    for p = 1:num_points
        % 使用对应颜色绘制
        plot(seq_list, squeeze(all_lamda(b,p,:)), ...
            '-o', ...
            'Color', my_colors{p}, ...
            'MarkerFaceColor', my_colors{p}, ...
            'MarkerSize', 5, ...
            'LineWidth', 1.5);
    end
    
    % 样式美化
    title(sprintf('%.1f T', B0_list(b)), 'FontSize', 12, 'FontWeight', 'bold');
    if b == 1
        ylabel({'\lambda^*';'(Sensitivity Var.)'}, 'FontSize', 11, 'FontWeight', 'bold');
    end
    xlabel(''); % 第一行不显示 X 轴标签，减少冗余
    xticklabels({}); 
    ylim(ylim_lamda);
    xlim([0 max(seq_list)+100]);
    grid on; box off;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
end

%% --- 3.3 绘制第二行：Alpha (Noise Variability) ---
for b = 1:numel(B0_list)
    ax = nexttile(t); hold on;
    
    for p = 1:num_points
        plot(seq_list, squeeze(all_alpha(b,p,:)), ...
            '-o', ...
            'Color', my_colors{p}, ...
            'MarkerFaceColor', my_colors{p}, ...
            'MarkerSize', 5, ...
            'LineWidth', 1.5);
    end
    
    % 在 7T 的 Alpha 图中 (最后一列)，标记截断点 (如果需要强调)
    if B0_list(b) == 7
        % 绘制一条虚线表示 N=2500 (或您讨论的 N=2000)
        xline(2500, '--k', 'Truncation', 'LabelVerticalAlignment','bottom');
    end

    if b == 1
        ylabel({'\alpha^*';'(Noise Var.)'}, 'FontSize', 11, 'FontWeight', 'bold');
    end
    xlabel('Number of Basis Modes');
    ylim(ylim_alpha);
    xlim([0 max(seq_list)+100]);
    grid on; box off;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
end

%% --- 3.4 添加全局图例 ---
% 使用 dummy plot 生成图例，放在图表底部
h = zeros(3, 1);
for p = 1:num_points
    h(p) = plot(nan, nan, '-o', 'Color', my_colors{p}, 'MarkerFaceColor', my_colors{p}, 'LineWidth', 1.5);
end
lg = legend(h, point_names, 'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south'; % 将图例置于整个布局的底部
lg.FontSize = 12;

% 总标题 (可选)
title(t, 'Convergence of Intrinsic Temporal Parameters (\lambda^*, \alpha^*)', ...
    'FontSize', 14, 'FontWeight', 'bold');