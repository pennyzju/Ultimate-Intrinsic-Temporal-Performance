clear; clc;

%% ================= 1. 参数与配置 =================
% 1.1 选择场强 (1=7T, 2=3T, 3=1.5T)
target_field_idx = 1; 

% 1.2 空间点坐标
points = [
    48 56 80;   % Center
    48 67 80;   % Intermediate
    48 75 80    % Edge
];
pointNames = {'Center Location', 'Intermediate Location', 'Edge Location'};
pointColors = [
    0.85, 0.33, 0.10; % Orange-Red
    0.93, 0.69, 0.13; % Yellow-Gold
    0.00, 0.45, 0.74  % Blue
];
num_points = size(points, 1);

% 1.3 场强配置
configs(1).name = '7 T';
configs(1).path = 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain'; 
configs(2).name = '3 T';
configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain';
configs(3).name = '1.5 T';
configs(3).path = 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain';

target_seq = 2500;
flag = 'left';

% === 修改点：大幅增加 K 的上限 ===
% 之前是 logspace(-1, 5, ...)，现在改为到 9 (即 10的9次方)
% 如果 lambda 很小 (比如 1e-4)，这就需要 snr 达到 1e4 才能看到效果
k_sim_vector = logspace(-1, 5, 100); 

%% ================= 2. 数据加载 =================
cfg = configs(target_field_idx);
fprintf('正在读取基准数据 (Field: %s, Seq: %d)...\n', cfg.name, target_seq);

snr_dir  = cfg.path; 
base_dir = fileparts(snr_dir);
loss_dir = fullfile(base_dir, 'loss_plain');
b1_dir   = fullfile(base_dir, 'B1_map');

base_alpha = zeros(num_points, 1);
base_lamda = zeros(num_points, 1);
base_snr   = zeros(num_points, 1); 

% --- 读取 Alpha ---
f_alpha = fullfile(loss_dir, sprintf('alpha_wo_%s_%d.mat', flag, target_seq));
if exist(f_alpha, 'file')
    tmp = load(f_alpha); 
    if isfield(tmp, 'alpha'), val = tmp.alpha; else, val = tmp.lamda; end
    val = abs(val);
    for p=1:num_points
        base_alpha(p) = val(points(p,1), points(p,2), points(p,3));
    end
end

% --- 读取 Lambda ---
f_lamda = fullfile(b1_dir, sprintf('lamda_%s_%d.mat', flag, target_seq));
if exist(f_lamda, 'file')
    tmp = load(f_lamda);
    if isfield(tmp, 'lamda'), val = tmp.lamda; else, val = tmp.alpha; end
    val = abs(val);
    for p=1:num_points
        base_lamda(p) = val(points(p,1), points(p,2), points(p,3));
    end
end

% --- 读取 Unit SNR ---
f_snr = fullfile(snr_dir, sprintf('USNR_left_-5_%d.mat', target_seq));
if ~exist(f_snr, 'file'), f_snr = fullfile(snr_dir, sprintf('B1_left_0_%d.mat', target_seq)); end
if exist(f_snr, 'file')
    tmp = load(f_snr); 
    vars = fields(tmp); val = tmp.(vars{1});
    for p=1:num_points
        base_snr(p) = val(points(p,1), points(p,2), points(p,3));
    end
end

%% ================= 3. 数据检查 (DEBUG) =================
fprintf('\n--- 检查点参数 ---\n');
for p = 1:num_points
    fprintf('%s: SNR_unit=%.4f, Alpha=%.4f, Lambda=%.4e\n', ...
        pointNames{p}, base_snr(p), base_alpha(p), base_lamda(p));
    if base_lamda(p) == 0
        warning('警告：Lambda 为 0，该曲线将永远是直线！请检查数据源。');
    end
end
fprintf('--------------------\n');

%% ================= 4. 优化绘图 (美化版) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]); % 窗口稍微宽一点
hold on;

% --- 1. 设置绘图范围 (关键步骤) ---
% 我们只关心 SNR > 0.1 的区域，因为那是饱和发生的地方
% 上限设为模拟的最大值
x_view_lims = [1e-1, max(k_sim_vector)*max(base_snr)]; 
% Y轴下限也设为 0.1，上限设为最大极限值的 2 倍，留出空间
max_limit = max(1./base_lamda(base_lamda>0));
y_view_lims = [1e-1, max_limit * 3]; 

% --- 2. 绘制参考线 ---
% 理想线只画在对角线上
plot(x_view_lims, x_view_lims, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
    'DisplayName', 'Ideal Linear (y=x)');

% --- 3. 循环绘制数据曲线 ---
for p = 1:num_points
    a = base_alpha(p);
    l = base_lamda(p);
    s = base_snr(p);
    
    % 计算
    input_snr_seq = k_sim_vector * s;
    denominator = sqrt(1 + a + (l .* input_snr_seq).^2);
    tsnr_seq = input_snr_seq ./ denominator;
    
    % --- 绘制主曲线 ---
    plot(input_snr_seq, tsnr_seq, '-', ...
        'Color', pointColors(p,:), ...
        'LineWidth', 3, ... % 线条加粗
        'DisplayName', pointNames{p});
    
    % --- 添加极限值标注 (Text Annotation) ---
    % 在曲线的最右端添加文字，说明极限 tSNR 是多少
    if l > 0
        limit_val = 1/l;
        
        % 画一条淡淡的水平渐近线
        yline(limit_val, ':', 'Color', pointColors(p,:), 'LineWidth', 1.2, 'Alpha', 0.6, 'HandleVisibility','off');
        
        % 在图表右侧添加数值文本
        text(x_view_lims(2)*0.6, limit_val*1.1, sprintf('Limit \\approx %.1f', limit_val), ...
            'Color', pointColors(p,:), ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'BackgroundColor', 'w', ... %以此遮挡背后的网格
            'Margin', 1);
    end
end

% --- 4. 坐标轴与装饰 ---
set(gca, 'XScale', 'log', 'YScale', 'log'); 
grid on; 
% 打开次级网格，但让它淡一点
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1); 

xlim(x_view_lims);
ylim(y_view_lims);

% 标签优化
xlabel('Total Input SNR (k \cdot SNR_{unit})', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('tSNR', 'FontSize', 13, 'FontWeight', 'bold');
title({['tSNR Saturation Curve (' cfg.name ')']; 'Zoomed in Saturation Region'}, 'FontSize', 15);

% 图例放在左上角，不仅不挡曲线，还利用了左上角的空白
legend('Location', 'northwest', 'FontSize', 12, 'Box', 'off');

set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Arial', 'TickDir', 'out');

fprintf('绘图完成。已裁剪左下角线性区域，聚焦饱和区。\n');

%% %% ================= 线性坐标绘图 (从 0 开始) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
hold on;

% --- 1. 重新定义 K 的范围 (线性生成，必须包含 0) ---
% 为了在线性图中看清弯曲过程，我们不需要模拟到 10^9 那么大
% 只需要模拟到 tSNR 基本饱和即可 (例如 3000 左右)
max_snr_view = 3000; 
% 反推需要的 k 值 (假设 base_snr 大约为 20-50 左右)
% k 从 0 开始，到大约 200
k_linear = [0, linspace(0.1, 200, 2000)]; 

% --- 2. 循环绘制 ---
for p = 1:num_points
    a = base_alpha(p);
    l = base_lamda(p);
    s = base_snr(p);
    
    % 计算 Input SNR (现在包含 0)
    input_snr_seq = k_linear * s;
    
    % 计算 tSNR
    % 当 input=0 时，分子为0，结果为0，符合逻辑
    denominator = sqrt(1 + a + (l .* input_snr_seq).^2);
    tsnr_seq = input_snr_seq ./ denominator;
    
    % 绘制
    plot(input_snr_seq, tsnr_seq, '-', ...
        'Color', pointColors(p,:), ...
        'LineWidth', 3, ...
        'DisplayName', pointNames{p});
    
    % --- 标注极限值 ---
    if l > 0
        limit_val = 1/l;
        yline(limit_val, ':', 'Color', pointColors(p,:), 'LineWidth', 1.5, 'HandleVisibility','off');
        
        % 在图表最右侧标注数值
        text(max_snr_view*0.9, limit_val - limit_val*0.05, sprintf('Max: %.1f', limit_val), ...
            'Color', pointColors(p,:), 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w');
    end
end

% --- 3. 绘制理想直线 y=x (作为参考) ---
% 在线性图中，y=x 会非常陡峭，我们只画一小段示意
plot([0, max_snr_view], [0, max_snr_view], '--', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Ideal (y=x)');


% --- 4. 关键设置：强制线性坐标并从 0 开始 ---
set(gca, 'XScale', 'linear', 'YScale', 'linear'); % 关闭对数！
xlim([0, max_snr_view]); 
ylim([0, max(1./base_lamda(base_lamda>0)) * 1.2]); % Y轴稍微留点空隙

grid on;
set(gca, 'GridAlpha', 0.4); 

xlabel('Total Input SNR (Linear Scale)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('tSNR', 'FontSize', 13, 'FontWeight', 'bold');
title({['tSNR Saturation (Linear View) - ' cfg.name]; 'Starting from SNR = 0'}, 'FontSize', 15);

legend('Location', 'southeast', 'FontSize', 12, 'Box', 'off'); % 图例改到右下角，因为左上角现在有曲线
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Arial');

fprintf('已切换为线性坐标，X轴起点为 0。\n');
