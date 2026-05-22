clear; clc; 
%% ================= 1. 参数与配置 =================
% 1.1 空间点坐标 [x, y, z]
points = [
    48 56 74;   % Center Location
    48 70 74;   % Intermediate Location (根据您的代码选的中间点)
    48 87 74    % Edge Location
];
pointNames = {'Center', 'Intermediate', 'Edge'};
num_points = size(points, 1);

% 1.2 场强配置 (结构体管理，防止顺序出错)
% 顺序：7T -> 3T -> 1.5T (对应 Input 2 的顺序)
configs(1).name = '7 T';
configs(1).path = 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain'; % 请确认路径
configs(1).scale = 49;  % SNR 乘数
configs(1).color = [0.85, 0.10, 0.10]; % 深红

configs(2).name = '3 T';
configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain';
configs(2).scale = 9;   % SNR 乘数
configs(2).color = [0.00, 0.60, 0.30]; % 深绿

configs(3).name = '1.5 T';
configs(3).path = 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain';
configs(3).scale = 2.25; % SNR 乘数
configs(3).color = [0.10, 0.30, 0.80]; % 深蓝

% 1.3 序列列表 (X轴)
seq_list = [20:20:100, 200, 400, 800, 1000, 1200, 1400,1500,1600,2000, 2200, 2500];
num_seqs = length(seq_list);
num_fields = length(configs);

% 文件名前缀/子文件夹定义
subdir_cv = '..'; % 注意：如果 CV 文件(alpha/lamda)不在 SNR_plain 里，需要回退或指定绝对路径
% 假设结构是：
% .../20251022seq/SNR_plain (存放 USNR)
% .../20251022seq/loss_plain (存放 alpha)
% .../20251022seq/B1_map     (存放 lamda)
flag = 'left'; 

%% ================= 2. 数据加载 =================
% 预分配矩阵 [场强, 空间点, 序列]
data_alpha = nan(num_fields, num_points, num_seqs);
data_lamda = nan(num_fields, num_points, num_seqs);
data_snr   = nan(num_fields, num_points, num_seqs);

fprintf('开始读取数据...\n');

for f = 1:num_fields
    % 推断路径
    snr_dir = configs(f).path; 
    % 假设 loss 和 B1 文件夹与 SNR_plain 同级
    base_dir = fileparts(snr_dir); % 回退一级
    loss_dir = fullfile(base_dir, 'loss_plain');
    b1_dir   = fullfile(base_dir, 'B1_map');
    
    fprintf('  Processing %s...\n', configs(f).name);
    
    for i = 1:num_seqs
        s_idx = seq_list(i);
        
        % --- 读取 Alpha ---
        f_alpha = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, s_idx));
        if exist(f_alpha, 'file')
            tmp = load(f_alpha); 
            % 兼容字段名可能是 alpha 或 lamda
            if isfield(tmp, 'alpha'), val = tmp.alpha; else, val = tmp.lamda; end
            val = abs(val);
            for p=1:num_points
                data_alpha(f, p, i) = val(points(p,1), points(p,2), points(p,3));
            end
        end
        
        % --- 读取 Lambda ---
        f_lamda = fullfile(b1_dir, sprintf('lamda_%s_%d.mat', flag, s_idx));
        if exist(f_lamda, 'file')
            tmp = load(f_lamda);
            if isfield(tmp, 'lamda'), val = tmp.lamda; else, val = tmp.alpha; end
            val = abs(val);
            for p=1:num_points
                data_lamda(f, p, i) = val(points(p,1), points(p,2), points(p,3));
            end
        end
        
        % --- 读取 SNR (并应用缩放) ---
        f_snr = fullfile(snr_dir, sprintf('USNR_left_0_%d.mat', s_idx));
        % 如果文件名是 B1_left... 请在这里修改
        if ~exist(f_snr, 'file')
             f_snr = fullfile(snr_dir, sprintf('B1_left_0_%d.mat', s_idx));
        end
        
        if exist(f_snr, 'file')
            tmp = load(f_snr); 
            vars = fields(tmp); val = tmp.(vars{1});
            for p=1:num_points
                raw_val = val(points(p,1), points(p,2), points(p,3));
                data_snr(f, p, i) = raw_val * configs(f).scale;
            end
        end
    end
end

% --- 对 SNR 取对数 ---
data_snr_log = log10(data_snr);

%% ================= 3. 计算统一纵坐标范围 =================
% 辅助函数：计算带边距的 limit
calc_lim = @(d) [min(d(:), [], 'omitnan') - range(d(:))*0.1, max(d(:), [], 'omitnan') + range(d(:))*0.1];

lim_alpha = calc_lim(data_alpha);
lim_lamda = calc_lim(data_lamda);
lim_snr   = calc_lim(data_snr_log);

% 防止全 NaN 或 range 为 0 的情况
if any(isnan(lim_alpha)), lim_alpha = [0 1]; end
if any(isnan(lim_lamda)), lim_lamda = [0 1]; end
if any(isnan(lim_snr)),   lim_snr   = [0 1]; end

%% ================= 4. 绘图 (3x3 布局) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.8, 0.8]);
% 3行3列，紧凑布局
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Row 1: SNR (Log10) ===
for p = 1:num_points
    nexttile; hold on;
    for f = 1:num_fields
        % 这里的 plot_data 已经是 log10 后的
        pl = plot(seq_list, squeeze(data_snr_log(f, p, :)), ...
            '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
        set(pl, 'MarkerFaceColor', pl.Color);
    end
    format_subplot(gca, lim_snr, seq_list);
    
    if p == 1, ylabel('Log_{10}(SNR)', 'FontWeight', 'bold'); end
    % 仅第一行显示列标题 (Location)
    title(pointNames{p}, 'FontSize', 12);
end


% === Row 2: Lambda ===
for p = 1:num_points
    nexttile; hold on;
    for f = 1:num_fields
        pl = plot(seq_list, squeeze(data_lamda(f, p, :)), ...
            '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
        set(pl, 'MarkerFaceColor', pl.Color);
    end
    format_subplot(gca, lim_lamda, seq_list);
    
    if p == 1, ylabel('\lambda (Sensitivity Var.)', 'FontWeight', 'bold'); end
end

% === Row 3:  Alpha===
for p = 1:num_points
    nexttile; hold on;
    for f = 1:num_fields
        pl = plot(seq_list, squeeze(data_alpha(f, p, :)), ...
            '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
        set(pl, 'MarkerFaceColor', pl.Color);
    end
    format_subplot(gca, lim_alpha, seq_list);
    
    % 仅第一列显示 Y 轴标签
    if p == 1, ylabel('\alpha (Noise Var.)', 'FontWeight', 'bold'); end
    
    xlabel('Number of Basis Modes'); % 只有最后一行显示 X Label
end

% === 共享图例 ===
% 创建虚拟图例对象
lgd = legend({configs.name}, 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south'; % 放在最底部

% 总标题
title(t, 'Convergence Comparison: SNR, \lambda, and \alpha across Field Strengths', ...
    'FontSize', 14, 'FontWeight', 'bold');


%% ================= 5. 辅助函数：统一格式化 =================
function format_subplot(ax, ylim_val, x_data)
    grid off;      % 无网格
    box off;       % 无边框
    xlim([0, max(x_data)+100]); % 统一 X 轴
    ylim(ylim_val);             % 统一 Y 轴
    
    set(ax, 'LineWidth', 1.2, 'FontSize', 14, 'TickDir', 'out', 'FontName', 'Arial');
    
    % 手动设置 Y 轴刻度，使其看起来不那么密
    % 如果是 SNR (log)，通常可以隔 0.5 或 1
    % 这里做一个简单的自适应刻度
    % y_ticks = linspace(ylim_val(1), ylim_val(2), 5);
    % set(ax, 'YTick', y_ticks);
end