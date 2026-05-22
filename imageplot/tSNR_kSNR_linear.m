clear; clc; close all;

%% ================= 1. 参数与配置 =================
% 定义三个感兴趣的点
points = [48 56 80; 48 70 80; 48 84 80];
pointNames = {'Center', 'Intermediate', 'Edge'};
%pointColors = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.00 0.45 0.74];
pointColors = [
    0.25 0.50 0.85;   % 柔和蓝
    0.20 0.65 0.35;   % 柔和绿
    0.80 0.30 0.30    % 柔和红
];
num_points = size(points, 1);

% 原始配置 (保持你的路径不变)
configs(1).name = '7 T';   configs(1).path = 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain'; 
configs(2).name = '3 T';   configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain';
configs(3).name = '1.5 T'; configs(3).path = 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain';

% 绘图顺序：从低场到高场 (1.5T -> 3T -> 7T)
% 对应 configs 数组的索引分别是 3, 2, 1
plot_order = [3, 2, 1]; 

target_seq = 2500;
flag = 'left';

%% ================= 2. 预读取数据以确定统一坐标轴 =================
fprintf('正在预读取所有数据以统一坐标范围...\n');
data_store = struct(); % 存储读取的数据
global_max_limit = 0;  % 用于记录所有情况下的最高 tSNR 极限

for i = 1:3
    cfg_idx = plot_order(i);
    cfg = configs(cfg_idx);
    
    snr_dir = cfg.path; 
    base_dir = fileparts(snr_dir);
    loss_dir = fullfile(base_dir, 'loss_plain');
    b1_dir   = fullfile(base_dir, 'B1_map');
    
    % 初始化数组
    curr_alpha = zeros(num_points, 1);
    curr_lamda = zeros(num_points, 1);
    
    % 读取 alpha
    f_alpha = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, target_seq));
    if exist(f_alpha,'file')
        tmp=load(f_alpha); 
        if isfield(tmp,'alpha'), v=tmp.alpha; else, v=tmp.lamda; end; v=abs(v);
        for p=1:num_points, curr_alpha(p)=v(points(p,1), points(p,2), points(p,3)); end
    else
        warning('文件未找到: %s', f_alpha);
    end
    
    % 读取 lambda
    f_lamda = fullfile(b1_dir, sprintf('lamda_%s_%d.mat', flag, target_seq));
    if exist(f_lamda,'file')
        tmp=load(f_lamda); 
        if isfield(tmp,'lamda'), v=tmp.lamda; else, v=tmp.alpha; end; v=abs(v);
        for p=1:num_points, curr_lamda(p)=v(points(p,1), points(p,2), points(p,3)); end
    else
        warning('文件未找到: %s', f_lamda);
    end
    
    % 存入结构体
    data_store(i).alpha = curr_alpha;
    data_store(i).lamda = curr_lamda;
    data_store(i).name  = cfg.name;
    
    % 更新全局最大值 (用于Y轴限制)
    % 理论极限 = 1 / lambda
    for p = 1:num_points
        if curr_lamda(p) > 1e-9 % 避免除以0
            limit_val = 1 / curr_lamda(p);
            if limit_val > global_max_limit
                global_max_limit = limit_val;
            end
        end
    end
end

% 设置防错默认值
if global_max_limit == 0 || isinf(global_max_limit)
    global_max_limit = 5000; 
end

% 定义统一的 X 轴范围 (通常设为最大极限的 3-5 倍以展示平台期)
x_limit_unified = global_max_limit * 4;
%x_limit_unified = 200;
if x_limit_unified > 60000, x_limit_unified = 60000; end % 封顶防止过大
if x_limit_unified < 500, x_limit_unified = 200; end

%% ================= 3. 绘图 (3个子图) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.3, 0.8, 0.45]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'tSNR Saturation Curves Across Field Strengths', 'FontSize', 16, 'FontWeight', 'bold');

% 生成 X 轴模拟数据
input_snr_sim = linspace(0, x_limit_unified, 2000);
y_limit_unified = global_max_limit * 1.1; % Y轴留 10% 顶部空间

legend_lines = []; % 用于存储图例句柄

for i = 1:3
    nexttile;
    hold on;
    
    curr_data = data_store(i);
    
    % 1. 绘制理想线 (y=x)
    h_ideal = plot([0, x_limit_unified], [0, x_limit_unified], '--', ...
        'Color', [0.8 0.8 0.8], 'LineWidth', 2.5, 'DisplayName', 'Ideal Linear');
    
    % 2. 绘制数据曲线
    for p = 1:num_points
        a = curr_data.alpha(p);
        l = curr_data.lamda(p);
        
        % tSNR 计算
        denominator = sqrt(1 + a + (l .* input_snr_sim).^2);
        tsnr_sim = input_snr_sim ./ denominator;
        
        % 绘制曲线
        h_line = plot(input_snr_sim, tsnr_sim, '-', ...
            'Color', pointColors(p,:), 'LineWidth', 3, 'DisplayName', pointNames{p});
        
        % 保存句柄用于生成统一图例 (只存最后一张图的即可)
        if i == 3
            legend_lines = [legend_lines; h_line];
        end
        
        % 3. 绘制极限虚线 (不带文字)
        if l > 0
            limit_val = 1/l;
            yline(limit_val, ':', 'Color', pointColors(p,:), 'LineWidth', 1.2, 'HandleVisibility','off');
        end
    end
    
    % 4. 设置坐标轴
    set(gca, 'XScale', 'linear', 'YScale', 'linear');
    xlim([0, x_limit_unified]);
    ylim([0, y_limit_unified]);
   % ylim([0,200])
    
   % grid on;
    set(gca, 'GridAlpha', 0.4, 'LineWidth', 1.2, 'FontSize', 12);
    
    xlabel('Total Input SNR', 'FontWeight', 'bold');
    if i == 1
        ylabel('tSNR', 'FontWeight', 'bold');
    end
    title(curr_data.name, 'FontSize', 14); % 显示 1.5T / 3T / 7T
end

%% ================= 4. 设置全局图例 =================
% 手动构建图例对象列表 (Ideal + 3个点)
all_legends = [h_ideal; legend_lines];
lgd = legend(all_legends, 'Orientation', 'horizontal', 'FontSize', 11);
lgd.Layout.Tile = 'south'; % 将图例放置在整个布局的最右侧

fprintf('绘图完成。\n统一 Y 轴范围: 0 - %.0f\n统一 X 轴范围: 0 - %.0f\n', y_limit_unified, x_limit_unified);