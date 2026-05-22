%% ================= 单独绘制 Alpha 比值 (Ratio) =================
% 确保之前的 data_alpha 存在
if ~exist('data_alpha', 'var')
    error('请先运行之前的代码以加载 data_alpha 数据！');
end

% 1. 计算比值
% data_alpha 维度: [场强(1=7T, 2=3T, 3=1.5T), 空间点, 序列]
ratio_7_3   = squeeze(data_alpha(1, :, :) ./ data_alpha(2, :, :)); % 7T / 3T
ratio_7_1p5 = squeeze(data_alpha(1, :, :) ./ data_alpha(3, :, :)); % 7T / 1.5T

% 2. 设置绘图参数
% 理论参考值 (假设 Alpha 与 B0 平方成正比)
ref_7_3   = (7/3)^2;    % ≈ 5.44
ref_7_1p5 = (7/1.5)^2;  % ≈ 21.78

ratio_configs(1).name  = 'Ratio: 7T / 3T';
ratio_configs(1).color = [0.85, 0.50, 0.10]; % 橙色
ratio_configs(1).ref   = ref_7_3;

ratio_configs(2).name  = 'Ratio: 7T / 1.5T';
ratio_configs(2).color = [0.60, 0.20, 0.80]; % 紫色
ratio_configs(2).ref   = ref_7_1p5;

% Y轴范围 (根据您的数据调整，建议 [0, 30])
ylim_ratio = [0, 1000]; 

%% 3. 开始绘图
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.3, 0.8, 0.45]); % 宽扁布局
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'normal');

for p = 1:num_points
    nexttile; hold on; box off;
    
    % --- A. 绘制实测数据曲线 ---
    % 7T / 3T
    plot(seq_list, ratio_7_3(p, :), '-s', ...
        'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', ratio_configs(1).color, 'MarkerFaceColor', ratio_configs(1).color, ...
        'DisplayName', ratio_configs(1).name);
    
    % 7T / 1.5T
    plot(seq_list, ratio_7_1p5(p, :), '-^', ...
        'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', ratio_configs(2).color, 'MarkerFaceColor', ratio_configs(2).color, ...
        'DisplayName', ratio_configs(2).name);
    
%     % --- B. 绘制理论参考线 (虚线) ---
%     yline(ref_7_3, '--', sprintf('Theory (%.1f)', ref_7_3), ...
%         'Color', ratio_configs(1).color, 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
%     
%     yline(ref_7_1p5, '--', sprintf('Theory (%.1f)', ref_7_1p5), ...
%         'Color', ratio_configs(2).color, 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
%     
    % --- C. 格式化坐标轴 ---
    grid on; grid minor;
    xlim([0, max(seq_list)+100]);
    ylim(ylim_ratio);
    
    xlabel('Number of Basis Modes', 'FontSize', 11, 'FontWeight', 'bold');
    title(pointNames{p}, 'FontSize', 13, 'FontWeight', 'bold');
    
    % 仅第一列显示 Y 轴标签
    if p == 1
        ylabel('Ratio of \alpha (Noise Var.)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    set(gca, 'FontSize', 10, 'LineWidth', 1.2, 'TickDir', 'out', 'FontName', 'Arial');
end

% --- 4. 添加共享图例 ---
% 只取前两个句柄做图例 (忽略虚线)
h_legend = findobj(gca, 'Type', 'Line'); 
% findobj 顺序通常是倒序的，取最后绘制的两个(即实线)
lgd = legend([h_legend(end), h_legend(end-1)], ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
lgd.Layout.Tile = 'north'; % 图例放在顶部

title(t, 'Field Strength Dependence of Intrinsic Noise Variability (\alpha)', ...
    'FontSize', 15, 'FontWeight', 'bold');