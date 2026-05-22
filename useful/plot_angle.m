%% 可视化：固定 seq 和 flag，观察不同 angle 的影响
loss_dir = 'H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain';    % 替换为你的数据存放路径
flag = 'left';      % 指定 flag
target_seq = 100;   % 指定 seq
slice_idx = 48;     % 指定层面

angles = -5:5;      % 角度范围：-5, -4, ..., 0, ..., 5
num_plots = numel(angles);

% 预分配空间
plot_data = cell(num_plots, 1);
c_min = inf;
c_max = -inf;

fprintf('Loading files for seq %d, flag %s...\n', target_seq, flag);

%% 1. 数据加载与范围计算
for i = 1:num_plots
    curr_angle = angles(i);
    
    % 根据文件名规则：loss_flag_angle_seq_correction
    % 注意：根据你之前代码的逻辑，angle=0 时可能是 baseline 文件，
    % 这里假设所有角度（包括0）都有对应的 correction 文件。
    fname = sprintf('loss_%s_%d_%d_correction.mat', flag, curr_angle, target_seq);
    path = fullfile(loss_dir, fname);
    
    if exist(path, 'file')
        S = load(path);
        % 提取切面 (48, :, :)
        % 如果变量名是 loss_map_correction，请根据实际修改
        slice = squeeze(S.loss_map_correction(slice_idx, :, :));
        plot_data{i} = slice;
        
        % 更新全局颜色刻度范围
        valid_vals = slice(isfinite(slice));
        if ~isempty(valid_vals)
            c_min = min(c_min, min(valid_vals(:)));
            c_max = max(c_max, max(valid_vals(:)));
        end
    else
        warning('未找到文件: %s', fname);
        plot_data{i} = [];
    end
end

%% 2. 绘图展示
figure('Color', 'w', 'Position', [50, 50, 1500, 600]);
% 使用 tiledlayout 布局 (2行6列 可容纳11-12个子图)
t = tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:num_plots
    nexttile;
    if isempty(plot_data{i})
        axis off; title(sprintf('Angle: %d (Missing)', angles(i)));
        continue;
    end
    
    imagesc(plot_data{i});
    clim([c_min, c_max]); % 统一颜色刻度
    colormap('jet');
    axis image off;
    title(sprintf('Angle: %d', angles(i)), 'FontSize', 10);
end

% 添加全局 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Correction Loss Value';

% 总标题
title(t, sprintf('Correction Maps (Flag: %s, Seq: %d, Slice: %d)', ...
    flag, target_seq, slice_idx), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Visualization complete.\n');