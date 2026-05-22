%% 形状位置对比脚本 (带基准边缘显示)
loss_dir = 'H:\UISNR\20240801_UISNR_output\20251022seq\loss_plain'; 
flag = 'left';
target_seq = 400;
slice_idx = 48;
angles = -5:5;

% --- 1. 处理基准数据 ---
baseline_name = sprintf('loss_%s_0_%d.mat', flag, target_seq);
S0 = load(fullfile(loss_dir, baseline_name));
loss0 = squeeze(S0.loss_map(slice_idx, :, :));

% 只要有数值就设为 1
mask0 = loss0 ~= 0; 

%% 2. 循环绘制 11 幅位置差异图
figure('Color', 'w', 'Position', [50, 50, 1600, 700]);
t = tiledlayout(2, 6, 'TileSpacing', 'none', 'Padding', 'compact');

for i = 1:numel(angles)
    curr_angle = angles(i);
    ax = nexttile;
    
    corr_name = sprintf('loss_%s_%d_%d_correction.mat', flag, curr_angle, target_seq);
    if exist(fullfile(loss_dir, corr_name), 'file')
        Sc = load(fullfile(loss_dir, corr_name));
        loss_curr_raw = squeeze(Sc.loss_map_correction(slice_idx, :, :));
        
        % 当前图二值化
        mask_curr = loss_curr_raw ~= 0;
        
        % 计算异或差异（不重合的地方为 1）
        location_diff = xor(mask0, mask_curr);
        
        % --- 绘图逻辑 ---
        % 1. 绘制差异背景：0显示深蓝，1显示红色
        imagesc(location_diff); 
        hold on;
        
        % 2. 绘制基准边缘：使用 contour 提取 mask0 的 0.5 分界线
        % 'g' 代表绿色边缘，你可以根据喜好换成 'w' (白色) 或 'y' (黄色)
        [~, hContour] = contour(mask0, [0.5 0.5], 'g', 'LineWidth', 1.2);
        
        axis image off;
        colormap(ax, [0 0 0.3; 1 0 0]); % 当前坐标系的色表：背景暗蓝，差异红
        
        % 计算偏移百分比
        base_area = sum(mask0(:));
        diff_ratio = (sum(location_diff(:)) / base_area) * 100;
        title(sprintf('Angle: %d (Diff: %.1f%%)', curr_angle, diff_ratio));
        hold off;
    else
        axis off;
    end
end

title(t, 'Location Mismatch (Red=Offset, Green=Baseline Edge)', 'FontSize', 14);