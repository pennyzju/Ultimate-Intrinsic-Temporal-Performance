%% CV 收敛对比曲线 (带下降点自动标注与双轴优化)
% ... (前面的数据加载代码保持不变) ...

%% 绘制
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% --- 左侧轴：绘制 Alpha (波动较大的 Loss CV) ---
yyaxis left
p1 = plot(seq_list, mean_alpha, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', '#0072BD');
ylabel('\alpha (Loss CV)', 'Color', '#0072BD', 'FontSize', 12, 'FontWeight', 'bold');
ax = gca; ax.YColor = '#0072BD'; % 轴颜色同步

% 自动寻找并标注 Alpha 的下降位置 (局部极小值)
[min_val, min_idx] = findpeaks(-mean_alpha); % 找负值的峰值即找原值的谷值
if ~isempty(min_idx)
    for k = 1:numel(min_idx)
        idx = min_idx(k);
        text(seq_list(idx), mean_alpha(idx), sprintf('  Dip: %d', seq_list(idx)), ...
            'Color', '#0072BD', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');
        plot(seq_list(idx), mean_alpha(idx), 'ro', 'MarkerSize', 12, 'LineWidth', 2); % 红圈标注
    end
end

% --- 右侧轴：绘制 Lamda (平稳下降的 B1 CV) ---
yyaxis right
p2 = plot(seq_list, mean_lamda, '-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', '#D95319');
ylabel('\lambda (B1 CV)', 'Color', '#D95319', 'FontSize', 12, 'FontWeight', 'bold');
ax.YColor = '#D95319'; % 轴颜色同步

% 通用修饰
grid on;
xlabel('Sequence Index (seq)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Convergence Comparison (Slice ', num2str(slice_idx), ')'], 'FontSize', 14);
legend([p1, p2], {'\alpha (Loss CV) - Left Axis', '\lambda (B1 CV) - Right Axis'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');

% 调整坐标范围让曲线分得更开，方便观察
yyaxis left; ylim([min(mean_alpha)*0.8, max(mean_alpha)*1.2]);
yyaxis right; ylim([min(mean_lamda)*0.8, max(mean_lamda)*1.2]);

hold off;