clear; clc; close all;

%% ================= 1. 基础配置 =================
a = 2; 
points = [48 56 80; 48 70 80; 48 83 80];
pointNames = {'Center', 'Intermediate', 'Edge'};
num_points = size(points, 1);

configs(1) = struct('name', '7 T',   'scale', 49,   'color', [0.85, 0.10, 0.10], 'path', 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain');
configs(2) = struct('name', '3 T',   'scale', 9,    'color', [0.00, 0.60, 0.30], 'path', 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain');
configs(3) = struct('name', '1.5 T', 'scale', 2.25, 'color', [0.10, 0.30, 0.80], 'path', 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain');

K_FACTOR = 10^(a+4); 
seq_list = [20:20:100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200, 2500];
num_seqs = length(seq_list);
num_fields = length(configs);
flag = 'left';

%% ================= 2. 数据加载与统一校准 =================
% 预分配干净的矩阵
data_snr    = nan(num_fields, num_points, num_seqs);
data_alpha  = nan(num_fields, num_points, num_seqs);
data_lambda = nan(num_fields, num_points, num_seqs);

for f = 1:num_fields
    snr_dir = configs(f).path;
    base_dir = fileparts(snr_dir);
    loss_dir = fullfile(base_dir, 'loss_plain');
    b1_dir   = fullfile(base_dir, 'B1_map');
    
    fprintf('Loading %s data...\n', configs(f).name);
    
    for i = 1:num_seqs
        s_idx = seq_list(i);
        
        % 1. 加载 SNR 并立即应用所有缩放 (关键优化：此处一次算清)
        f_snr = fullfile(snr_dir, sprintf('USNR_left_0_%d.mat', s_idx));
        if ~exist(f_snr, 'file'), f_snr = fullfile(snr_dir, sprintf('B1_left_0_%d.mat', s_idx)); end
        
        if exist(f_snr, 'file')
            tmp = load(f_snr);
            vars = fields(tmp); 
            val_map = tmp.(vars{1});
            for p = 1:num_points
                % 最终物理 SNR = 原始值 * scale * K
                data_snr(f, p, i) = val_map(points(p,1), points(p,2), points(p,3)) * configs(f).scale * K_FACTOR;
            end
        end

        % 2. 加载 Alpha
        f_alpha = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, s_idx));
        if exist(f_alpha, 'file')
            tmp = load(f_alpha);
            val_map = if_field_else(tmp, 'alpha', 'lamda');
            for p = 1:num_points
                data_alpha(f, p, i) = abs(val_map(points(p,1), points(p,2), points(p,3)));
            end
        end

        % 3. 加载 Lambda
        f_lamda = fullfile(b1_dir, sprintf('lamda_%s_%d.mat', flag, s_idx));
        if exist(f_lamda, 'file')
            tmp = load(f_lamda);
            val_map = if_field_else(tmp, 'lamda', 'alpha');
            for p = 1:num_points
                data_lambda(f, p, i) = abs(val_map(points(p,1), points(p,2), points(p,3)));
            end
        end
    end
end

%% ================= 3. 核心指标计算 =================
% 此时所有输入量 (data_snr, data_alpha, data_lambda) 已经是最终量级
data_tsnr   = data_snr ./ sqrt(1 + data_alpha + (data_lambda .* data_snr).^2);
data_sratio = data_tsnr ./ data_snr;
data_snr_log = log10(data_snr);

%% ================= 4. 绘图 =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.8]);
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘图配置表：{数据源, Y轴标签, Y轴范围}
plot_configs = {
    data_snr,    'SNR',               calc_lim(data_snr);
    data_tsnr,   'tSNR',              calc_lim(data_tsnr);
    data_sratio, 'S_{ratio} (t/s)',   [0, 1.1]
};

for row = 1:3
    cur_data = plot_configs{row, 1};
    for p = 1:num_points
        ax = nexttile; hold on;
        for f = 1:num_fields
            plot(seq_list, squeeze(cur_data(f, p, :)), '-o', ...
                'LineWidth', 1.2, 'MarkerSize', 4, ...
                'Color', configs(f).color, 'MarkerFaceColor', configs(f).color);
        end
        format_subplot(ax, plot_configs{row, 3}, seq_list);
        if p == 1, ylabel(plot_configs{row, 2}, 'FontWeight', 'bold'); end
        if row == 1, title(pointNames{p}); end
        if row == 3, xlabel('Number of Basis Vectors'); end
    end
end
lgd = legend({configs.name}, 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';

%% ================= 4. 绘图 (对数展示与手动范围控制) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.6]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 【手动指定 Y 轴范围】在这里修改 [min, max]
% 如果设为 []，则会自动调用 calc_lim 函数计算自适应范围
y_lim_snr  = [0.6,4];  % 例如：log10(SNR) 的范围
y_lim_tsnr = [0.6,4];    % 例如：log10(tSNR) 的范围，1.3 对应 log10(20)

% 绘图配置表：{数据源(取log), Y轴标签, Y轴范围}
plot_configs = {
    log10(data_snr),  'log_{10}(SNR)', y_lim_snr;
    log10(data_tsnr), 'log_{10}(tSNR)', y_lim_tsnr
};

for row = 1:2
    cur_data = plot_configs{row, 1};
    y_label_str = plot_configs{row, 2};
    manual_lim = plot_configs{row, 3};
    
    % 如果未手动指定(为空)，则自动计算
    if isempty(manual_lim)
        y_lims = calc_lim(cur_data);
    else
        y_lims = manual_lim;
    end

    for p = 1:num_points
        ax = nexttile; hold on;
        for f = 1:num_fields
            % 提取当前场强、当前位置随序列变化的向量
            y_vals = squeeze(cur_data(f, p, :));
            
            % 绘图
            p_obj = plot(seq_list, y_vals, '-o', ...
                'LineWidth', 1.3, 'MarkerSize', 5, ...
                'Color', configs(f).color, ...
                'MarkerFaceColor', configs(f).color);
        end
        
        % 格式化子图
        format_subplot(ax, y_lims, seq_list);
        
        % 标签处理
        if p == 1, ylabel(y_label_str, 'FontWeight', 'bold', 'FontSize', 20); end
        if row == 1, title(pointNames{p}, 'FontSize', 13); end
        if row == 2, xlabel('Number of Basis Vectors', 'FontSize', 11); end
    end
end

% 共享图例
lgd = legend({configs.name}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 18);
lgd.Layout.Tile = 'south';

%% ================= 5. 辅助函数 (保持不变) =================
function lims = calc_lim(d)
    v = d(:); v = v(~isnan(v) & ~isinf(v));
    if isempty(v), lims = [0 1]; return; end
    margin = (max(v) - min(v)) * 0.1;
    if margin == 0, margin = 0.5; end
    lims = [min(v) - margin, max(v) + margin];
end

function format_subplot(ax, ylim_val, x_data)
    grid off; box off; 
    set(ax, 'XGrid', 'off', 'YGrid', 'off', 'GridLineStyle', ':', 'GridAlpha', 0.5);
    xlim([0, max(x_data)*1.05]); 
    ylim(ylim_val);
    set(ax, 'LineWidth', 1.1, 'FontSize', 14, 'TickDir', 'in', 'FontName', 'Arial');
end

function val = if_field_else(s, f1, f2)
    if isfield(s, f1), val = s.(f1); else, val = s.(f2); end
end




