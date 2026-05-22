clear; clc; close all;
%暂时不需要合并

%% ================= 1. 参数与配置 =================
% 1.1 空间点坐标 [x, y, z]
points = [
    48 56 80;   % Center Location
    48 70 80;   % Intermediate Location
    48 74 80    % Edge Location
];
pointNames = {'Center Location', 'Intermediate Location', 'Edge Location'};
num_points = size(points, 1);

% 1.2 场强配置 (顺序：7T -> 3T -> 1.5T)
configs(1).name = '7 T';
configs(1).path = 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain'; 
configs(1).scale = 49;  
configs(1).color = [0.85, 0.10, 0.10]; % 深红

configs(2).name = '3 T';
configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain';
configs(2).scale = 9;   
configs(2).color = [0.00, 0.60, 0.30]; % 深绿

configs(3).name = '1.5 T';
configs(3).path = 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain';
configs(3).scale = 2.25; 
configs(3).color = [0.10, 0.30, 0.80]; % 深蓝

% 1.3 序列列表 (X轴)
seq_list = [ 20:20:100, 200, 400, 800, 1000, 1200, 1400,1500,1600, 2000, 2200, 2500];
num_seqs = length(seq_list);
num_fields = length(configs);

% 文件名前缀/子文件夹定义
subdir_cv = '..'; 
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
    base_dir = fileparts(snr_dir); 
    loss_dir = fullfile(base_dir, 'loss_plain');
    b1_dir   = fullfile(base_dir, 'B1_map');
    
    fprintf('  Processing %s...\n', configs(f).name);
    
    for i = 1:num_seqs
        s_idx = seq_list(i);
        
        % --- 读取 Alpha ---
        f_alpha = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, s_idx));
        if exist(f_alpha, 'file')
            tmp = load(f_alpha); 
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

%% ================= 3. 计算比值 (Ratio) =================
% 假设 configs(1)=7T, configs(2)=3T, configs(3)=1.5T
% 我们计算 7T/3T 和 7T/1.5T
% 结果矩阵维度: [比值类型(2), 空间点, 序列]
ratio_data = nan(2, num_points, num_seqs);

% Ratio 1: 7T / 3T
ratio_data(1, :, :) = data_alpha(1, :, :) ./ data_alpha(2, :, :);

% Ratio 2: 7T / 1.5T
ratio_data(2, :, :) = data_alpha(2, :, :) ./ data_alpha(3, :, :);

% 定义比值的颜色和名称
ratio_configs(1).name = '7T / 3T';
ratio_configs(1).color = [0.85, 0.50, 0.10]; % 橙色
ratio_configs(2).name = '3T / 1.5T';
ratio_configs(2).color = [0.50, 0.00, 0.50]; % 紫色

%% ================= 4. 计算统一纵坐标范围 =================
calc_lim = @(d) [min(d(:), [], 'omitnan') - range(d(:))*0.1, max(d(:), [], 'omitnan') + range(d(:))*0.1];

lim_alpha = calc_lim(data_alpha);
lim_lamda = calc_lim(data_lamda);
lim_snr   = calc_lim(data_snr_log);
%lim_ratio = calc_lim(ratio_data); % Ratio 的统一范围
lim_ratio =[0,3];
% 安全检查
if any(isnan(lim_alpha)), lim_alpha = [0 1]; end
if any(isnan(lim_lamda)), lim_lamda = [0 1]; end
if any(isnan(lim_snr)),   lim_snr   = [0 1]; end
if any(isnan(lim_ratio)), lim_ratio = [0 1]; end

%% ================= 5. 绘图 (4x3 布局) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.6, 0.9]);
% 改为 4 行 3 列
t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% === Row 1: SNR (Log10) ===
for p = 1:num_points
    nexttile; hold on;
    for f = 1:num_fields
        pl = plot(seq_list, squeeze(data_snr_log(f, p, :)), ...
            '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
        set(pl, 'MarkerFaceColor', pl.Color);
    end
    format_subplot(gca, lim_snr, seq_list);
    if p == 1, ylabel('Log_{10}(SNR)', 'FontWeight', 'bold'); end
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
    if p == 1, ylabel('\lambda', 'FontWeight', 'bold'); end
end

% === Row 3: Alpha ===
for p = 1:num_points
    nexttile; hold on;
    for f = 1:num_fields
        pl = plot(seq_list, squeeze(data_alpha(f, p, :)), ...
            '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
        set(pl, 'MarkerFaceColor', pl.Color);
    end
    format_subplot(gca, lim_alpha, seq_list);
    if p == 1, ylabel('\alpha', 'FontWeight', 'bold'); end
end

% === Row 4: Alpha Ratio (New!) ===
for p = 1:num_points
    nexttile; hold on;
    % 绘制 7T/3T
    pl1 = plot(seq_list, squeeze(ratio_data(1, p, :)), ...
        '-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', ratio_configs(1).color);
    set(pl1, 'MarkerFaceColor', pl1.Color);
    
    % 绘制 7T/1.5T
    pl2 = plot(seq_list, squeeze(ratio_data(2, p, :)), ...
        '-^', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', ratio_configs(2).color);
    set(pl2, 'MarkerFaceColor', pl2.Color);
    
    format_subplot(gca, lim_ratio, seq_list);
    xlabel('Number of Basis Modes'); 
    if p == 1, ylabel('\alpha Ratio', 'FontWeight', 'bold'); end
    
    % 可选：加一条 y=1 的参考线
    yline(1, '--k', 'LineWidth', 0.8, 'Alpha', 0.3);
end

% === 图例设置 ===
% 由于上面 3 行和第 4 行的图例内容不一样，我们可以用一个 Trick
% 在底部放两个 Legend，或者合并在一起

% 方案：只在底部放总图例
% 构造 dummy lines 用于图例
h = [];
% 1. 场强图例
for f = 1:num_fields
    h(end+1) = plot(nan, nan, '-o', 'Color', configs(f).color, 'LineWidth', 2, 'DisplayName', configs(f).name);
end
% 2. Ratio 图例
for r = 1:2
    h(end+1) = plot(nan, nan, '-s', 'Color', ratio_configs(r).color, 'LineWidth', 2, 'DisplayName', ratio_configs(r).name);
end

lgd = legend(h, 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';

title(t, 'Convergence & Field Ratio Analysis', 'FontSize', 14, 'FontWeight', 'bold');


%% ================= 6. 辅助函数 =================
function format_subplot(ax, ylim_val, x_data)
    %grid on;       % 开启网格方便看比值
    %grid minor;
    box off;
    xlim([0, max(x_data)+100]); 
    % 为防止 Y 轴范围过大或过小，做个简单判断
    if ylim_val(2) > ylim_val(1)
        ylim(ylim_val); 
    end
    set(ax, 'LineWidth', 1.0, 'FontSize', 9, 'TickDir', 'out', 'FontName', 'Arial');
end