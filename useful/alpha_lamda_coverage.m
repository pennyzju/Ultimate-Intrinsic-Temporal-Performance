%% CV 收敛对比曲线 (Alpha vs Lamda)
B0 =3;
subdir = '20251022seq';
loss_dir = fullfile( B0_folderpath(B0),subdir,'loss_plain');
B1_dir = fullfile( B0_folderpath(B0),subdir,'B1_map');
flag = 'left';
seq_list = [100, 200, 400, 800, 1000, 1200, 1500,2000,2200];
slice_idx = 48;

mean_alpha = zeros(size(seq_list));
mean_lamda = zeros(size(seq_list));

for i = 1:numel(seq_list)
    curr_seq = seq_list(i);
    
    % --- 1. 加载 Alpha 数据 ---
    alpha_file = fullfile(loss_dir, sprintf('alpha_%s_%d.mat', flag, curr_seq));
    if exist(alpha_file, 'file')
        a_data = load(alpha_file);
        % 兼容不同变量名，并取绝对值处理复数
        if isfield(a_data, 'alpha'), a_val = a_data.alpha; else, a_val = a_data.lamda; end
        temp_a = abs(squeeze(a_val(slice_idx, :, :)));
        mean_alpha(i) = mean(temp_a(temp_a > 0 & isfinite(temp_a)), 'all');
    else
        mean_alpha(i) = NaN;
    end
    
    % --- 2. 加载 Lamda 数据 ---
    lamda_file = fullfile(B1_dir, sprintf('lamda_%s_%d.mat', flag, curr_seq));
    if exist(lamda_file, 'file')
        l_data = load(lamda_file);
        if isfield(l_data, 'lamda'), l_val = l_data.lamda; else, l_val = l_data.alpha; end
        temp_l = abs(squeeze(l_val(slice_idx, :, :)));
        mean_lamda(i) = mean(temp_l(temp_l > 0 & isfinite(temp_l)), 'all');
    else
        mean_lamda(i) = NaN;
    end
end

%% 3. 绘图
figure('Color', 'w', 'Position', [200, 200, 800, 500]);
hold on;

% 绘制 Alpha 曲线 (蓝色)
p1 = plot(seq_list, mean_alpha, '-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', '#0072BD', 'MarkerFaceColor', '#0072BD');

% 绘制 Lamda 曲线 (橙色)
p2 = plot(seq_list, mean_lamda, '-s', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Color', '#D95319', 'MarkerFaceColor', '#D95319');

% 修饰图表
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1);
xlabel('Sequence Index (seq)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coefficient of Variation (CV)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Convergence Comparison: \\alpha vs \\lambda (Slice %d)', slice_idx), 'FontSize', 14);

% 添加图例
legend([p1, p2], {'\alpha (Loss CV)', '\lambda (B1 CV)'}, 'Location', 'northeast');

% 如果量级差异巨大，取消下面两行的注释使用对数坐标
% set(gca, 'YScale', 'log');

hold off;