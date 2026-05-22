clear; clc; close all;

%% ================= 1. 基础配置 =================
a = 9; 
% 选取的点（增加至 6 个，你可以根据实际坐标修改）
points = [48 56 80; 48 70 80; 48 83 80; 50 50 80; 40 40 80; 60 60 80];
pointNames = {'Center', 'Inter-1', 'Edge', 'Pos-4', 'Pos-5', 'Pos-6'};
num_points = size(points, 1);

% 场强配置 (横坐标数据)
field_strengths = [7, 3, 1.5]; % 用于横坐标刻度
configs(1) = struct('name', '7 T',   'field', 7,   'scale', 49,   'path', 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain');
configs(2) = struct('name', '3 T',   'field', 3,   'scale', 9,    'path', 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain');
configs(3) = struct('name', '1.5 T', 'field', 1.5, 'scale', 2.25, 'path', 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain');

K_FACTOR = 10^(a+4); 
flag = 'left';

% 选择一个特定的序列模式数进行对比 (例如：第 100 个模式)
target_seq = 100; 

%% ================= 2. 数据加载与计算 =================
% 结果存储：[场强数量 x 空间点数量]
res_snr  = nan(length(configs), num_points);
res_tsnr = nan(length(configs), num_points);

for f = 1:length(configs)
    fprintf('Processing %s for Field Strength analysis...\n', configs(f).name);
    
    % 路径构建
    base_dir = fileparts(configs(f).path);
    snr_path   = fullfile(configs(f).path, sprintf('USNR_left_0_%d.mat', target_seq));
    if ~exist(snr_path, 'file'), snr_path = fullfile(configs(f).path, sprintf('B1_left_0_%d.mat', target_seq)); end
    
    alpha_path = fullfile(base_dir, 'loss_plain', sprintf('alpha_wo_%s_%d.mat', flag, target_seq));
    lamda_path = fullfile(base_dir, 'B1_map', sprintf('lamda_%s_%d.mat', flag, target_seq));

    % 加载数据并计算 (一次性算清)
    if exist(snr_path, 'file') && exist(alpha_path, 'file') && exist(lamda_path, 'file')
        % 读取 SNR
        tmp = load(snr_path); vars = fields(tmp); 
        val_snr_map = tmp.(vars{1}) * configs(f).scale * K_FACTOR;
        
        % 读取 Alpha
        tmp = load(alpha_path); 
        if isfield(tmp, 'alpha'), val_alpha_map = abs(tmp.alpha); else, val_alpha_map = abs(tmp.lamda); end
        
        % 读取 Lambda
        tmp = load(lamda_path);
        if isfield(tmp, 'lamda'), val_lamda_map = abs(tmp.lamda); else, val_lamda_map = abs(tmp.alpha); end

        % 提取各点数值并计算指标
        for p = 1:num_points
            s = val_snr_map(points(p,1), points(p,2), points(p,3));
            a_val = val_alpha_map(points(p,1), points(p,2), points(p,3));
            l = val_lamda_map(points(p,1), points(p,2), points(p,3));
            
            t = s / sqrt(1 + a_val + (l * s)^2);
            
            res_snr(f, p) = s;
            res_tsnr(f, p) = t;
        end
    else
        warning('Data missing for %s at seq %d', configs(f).name, target_seq);
    end
end

%% ================= 3. 绘图 (横坐标为场强) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.5]);
tlo = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

% 手动范围设置 (如果不需要请设为 [])
y_lim_snr_log  = []; 
y_lim_tsnr_log = [];

% 绘图 1: log10(SNR)
nexttile; hold on;
color_map = lines(num_points); % 为不同点分配不同颜色
for p = 1:num_points
    plot(field_strengths, log10(res_snr(:, p)), '-s', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', color_map(p,:), 'Color', color_map(p,:), 'DisplayName', pointNames{p});
end
format_field_plot(gca, 'Field Strength (T)', 'log_{10}(SNR)', y_lim_snr_log, field_strengths);
grid on; legend('Location', 'northeastoutside');

% 绘图 2: log10(tSNR)
nexttile; hold on;
for p = 1:num_points
    plot(field_strengths, log10(res_tsnr(:, p)), '-d', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', color_map(p,:), 'Color', color_map(p,:), 'DisplayName', pointNames{p});
end
format_field_plot(gca, 'Field Strength (T)', 'log_{10}(tSNR)', y_lim_tsnr_log, field_strengths);
grid on;

title(tlo, sprintf('Comparison Across Field Strengths (at Modes = %d)', target_seq), 'FontSize', 14, 'FontWeight', 'bold');

%% ================= 辅助函数 =================
function format_field_plot(ax, x_lab, y_lab, y_lims, x_ticks)
    xlabel(x_lab, 'FontWeight', 'bold');
    ylabel(y_lab, 'FontWeight', 'bold');
    if ~isempty(y_lims), ylim(y_lims); end
    set(ax, 'XTick', sort(x_ticks), 'XDir', 'normal', 'FontSize', 11, 'LineWidth', 1.1);
    % 如果场强是从大到小排的，可以用 set(ax, 'XDir', 'reverse') 翻转坐标轴
end