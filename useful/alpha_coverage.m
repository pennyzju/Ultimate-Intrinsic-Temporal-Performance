%% 均值收敛曲线绘制
loss_dir = 'H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain'; % 文件夹路径
flag = 'left';
seq_list = [100, 200, 400, 600,800, 1000, 1200, 1400,1450,1500,1600,1700,1800, 2000];
slice_idx = 48; % 你关注的层面

mean_values = zeros(size(seq_list)); % 用于存储每个 seq 的均值

for i = 1:numel(seq_list)
    fname = sprintf('alpha_wo_%s_%d.mat', flag, seq_list(i));
    path = fullfile(loss_dir, fname);
    
    if exist(path, 'file')
        data = load(path);
        % 获取指定层面的 alpha
        temp_alpha = data.alpha(slice_idx, :, :);
        
        % 计算该层面的全局均值 (排除 NaN 和 Inf)
        valid_idx = isfinite(temp_alpha);
        if any(valid_idx(:))
            mean_values(i) = mean(temp_alpha(valid_idx));
        else
            mean_values(i) = NaN;
        end
    else
        mean_values(i) = NaN;
        warning('文件未找到: %s', fname);
    end
end

%% 绘图
figure('Color', 'w', 'Name', 'Convergence Curve');
hold on;

% 绘制折线图
plot(seq_list, mean_values, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', '#0072BD', 'Color', '#0072BD');

% 添加网格和标签
grid on;
xlabel('Sequence Index (seq)', 'FontSize', 12);
ylabel('Mean Alpha (CV)', 'FontSize', 12);
title(sprintf('Convergence of Alpha Mean (Slice %d)', slice_idx), 'FontSize', 14);

% 如果 seq 跨度很大，建议使用对数坐标
% set(gca, 'XScale', 'log'); 

% 可选：添加一条代表最终收敛值的参考线
line([seq_list(1) seq_list(end)], [mean_values(end) mean_values(end)], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
legend('Observed Mean', 'Final Level', 'Location', 'southeast');

hold off;