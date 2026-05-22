clear; clc; close all;

%% ================= 1. 参数与配置 =================
a =0; % 假设 a 的值，你可以根据实际需求修改
points = [48 56 80; 48 70 80; 48 83 80];
pointNames = {'Center', 'Intermediate', 'Edge'};
num_points = size(points, 1);

% 1.2 场强配置
configs(1).name = '7 T';
configs(1).path = 'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain';
configs(1).scale = 49;
configs(1).k =  10^(a+4); % 定义 k
configs(1).color = [0.85, 0.10, 0.10];

configs(2).name = '3 T';
configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain';
configs(2).scale = 9;
configs(2).k = 10^(a+4); % 定义 k
configs(2).color = [0.00, 0.60, 0.30];

configs(3).name = '1.5 T';
configs(3).path = 'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain';
configs(3).scale = 2.25; 
configs(3).k =  10^(a+4); % 定义 k
configs(3).color = [0.10, 0.30, 0.80];

seq_list = [20:20:100, 200, 400, 800, 1000, 1200, 1500,2000, 2200, 2500];
num_seqs = length(seq_list);
num_fields = length(configs);
flag = 'left'; 

%% ================= 2. 数据加载 (同前，略作精简说明) =================
data_alpha = nan(num_fields, num_points, num_seqs);
data_lamda = nan(num_fields, num_points, num_seqs);
data_snr_raw = nan(num_fields, num_points, num_seqs); % 存储未乘scale的原始SNR

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
                data_snr_raw(f, p, i) = raw_val * configs(f).scale* configs(f).k;
            end
        end
    end
end


%% ================= 3. 计算 tSNR 和 Sratio =================
data_tsnr = nan(num_fields, num_points, num_seqs);
data_sratio = nan(num_fields, num_points, num_seqs);

for f = 1:num_fields
    k = configs(f).k;
    for p = 1:num_points
        % 提取当前场强和位置的向量 [1 x num_seqs]
        snr_val = squeeze(data_snr_raw(f, p, :))'; 
        alpha_val = squeeze(data_alpha(f, p, :))';
        lamda_val = squeeze(data_lamda(f, p, :))';

        numerator =snr_val;
        denominator = sqrt(1 + alpha_val + (lamda_val .* numerator).^2);

        tsnr = numerator ./ denominator;
        sratio = tsnr ./ numerator;
        
        data_tsnr(f, p, :) = tsnr;
        data_sratio(f, p, :) = sratio;
    end
end

% 处理 SNR 的 Log 显示（使用原代码逻辑）
data_snr_final = nan(num_fields, num_points, num_seqs);
for f = 1:num_fields
    data_snr_final(f,:,:) = data_snr_raw(f,:,:) * configs(f).scale*configs(f).k;
end
data_snr_log = log10(data_snr_final);

%% ================= 4. 绘图 (扩展为 5x3 布局) =================
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.8, 0.9]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 准备绘图数据和标签的映射
plot_sets = {
    data_snr_final, 'SNR', calc_lim(data_snr_final);
    %data_snr_log, 'log_{10}(SNR)', calc_lim(data_snr_log);
    %data_lamda,   '\lambda (Sens. Var.)', calc_lim(data_lamda);
   %data_alpha,   '\alpha (Noise Var.)', calc_lim(data_alpha);
    data_tsnr,    'tSNR', calc_lim(data_tsnr);
    %data_sratio,  'S_{ratio}=tsnr/snr', [0, 1.1] % Ratio 通常在 0-1 之间
};

for row = 1:3
    current_data = plot_sets{row, 1};
    y_label_str = plot_sets{row, 2};
    y_lims = plot_sets{row, 3};
    
    for p = 1:num_points
        nexttile; hold on;
        for f = 1:num_fields
            pl = plot(seq_list, squeeze(current_data(f, p, :)), ...
                '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'Color', configs(f).color);
            set(pl, 'MarkerFaceColor', pl.Color);
        end
        format_subplot(gca, y_lims, seq_list);
        
        if p == 1, ylabel(y_label_str, 'FontWeight', 'bold'); end
        if row == 1, title(pointNames{p}, 'FontSize', 12); end
        if row == 2, xlabel('Number of Basis set'); end
    end
end

% 共享图例
lgd = legend({configs.name}, 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';

%% ================= 5. 辅助函数 =================
function lims = calc_lim(d)
    val_min = min(d(:), [], 'omitnan');
    val_max = max(d(:), [], 'omitnan');
    if isempty(val_min) || isnan(val_min), lims = [0 1]; return; end
    margin = (val_max - val_min) * 0.1;
    if margin == 0, margin = 0.1; end
    lims = [val_min - margin, val_max + margin];
end

function format_subplot(ax, ylim_val, x_data)
    grid off; box off;
    xlim([0, max(x_data)+100]);
    ylim(ylim_val);
    set(ax, 'LineWidth', 1.1, 'FontSize', 11, 'TickDir', 'out', 'FontName', 'Arial');
end