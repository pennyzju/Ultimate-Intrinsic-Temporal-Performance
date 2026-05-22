clear; clc; close all;

%% 1. 数据准备
% 每一行是一个体素坐标 [x y z]
points = [
    48 56 74;   % 中心
    %48 60 74;   % 右
    48 70 74; 
    %48 80 74;
    48 84 74;% 左
];

pointNames = {'Center Location', 'Intermediate Location', 'Edge Location'};

% 文件夹与参数
folders = {
    'H:\UISNR\20240801_UISNR_output\20251022seq\SNR_plain', ... % 7T
    'G:\20251201_3T_UISNR_output\20251022seq\SNR_plain',    ... % 3T
    'F:\20251201_1p5T_UISNR_output\20251022seq\SNR_plain'   ... % 1.5T
};
scaleFactors = [49, 9, 2.25];
folderNames = {'7 T', '3 T', '1.5 T'};
x_vals = [10:10:100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200, 2500];

% 颜色定义
sciColors = [
    0.85, 0.10, 0.10;  % 7T (深红)
    0.00, 0.60, 0.30;  % 3T (深绿)
    0.10, 0.30, 0.80   % 1.5T (深蓝)
];

%% 2. 读取所有数据并计算全局范围
num_files = length(x_vals);
num_points = size(points, 1);
num_folders = length(folders);
all_data = zeros(num_folders, num_points, num_files);

for f = 1:num_folders
    currentDir = folders{f};
    currentScale = scaleFactors(f);
    for i = 1:num_files
        % 确保文件名格式正确
        filename = fullfile(currentDir, sprintf('USNR_left_0_%d.mat', x_vals(i))); 
        if exist(filename, 'file')
            content = load(filename); vars = fields(content); data = content.(vars{1});
            for p = 1:num_points
                all_data(f, p, i) = data(points(p,1), points(p,2), points(p,3)) * currentScale;
            end
        else
            all_data(f, p, i) = NaN;
        end
    end
end
%% % === 关键步骤：计算全局统一的 Y 轴范围 ===

% 1. 取对数
log_data_all = log10(all_data);

% 2. 【核心修改】提取所有“有限值” (即剔除 NaN 和 +/-Inf)
% isfinite 会返回一个逻辑矩阵，只保留正常的数值
valid_vals = log_data_all(isfinite(log_data_all));

% 3. 安全计算 Min/Max
if isempty(valid_vals)
    % 如果数据全是 NaN 或 0，防止报错，给一个默认范围
    warning('数据全为无效值 (NaN/Inf)，使用默认坐标轴范围');
    unified_ylim = [0 1]; 
else
    global_min = min(valid_vals);
    global_max = max(valid_vals);

    % 4. 计算余量 (Padding)
    % 如果 max 和 min 相等 (直线)，强制给一个宽度
    if global_max == global_min
        padding = 1; 
    else
        padding = (global_max - global_min) * 0.05;
    end
    
    unified_ylim = [global_min - padding, global_max + padding];
end


%% 3. 绘图
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.5]); 
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:num_points
    nexttile; % 激活子图
    hold on; 
    
    % 【修改点1 & 2】: 
    % grid off (默认就是 off，这里确保一下)
    % box off (去掉右边和上边的黑框，只留坐标轴)
    grid off; 
    box off; 
    
    for f = 1:num_folders
        y_data = squeeze(all_data(f, p, :));
        plot_y = log10(y_data); % 取对数
        
        plot(x_vals, plot_y, ...
            'LineStyle', '-', 'LineWidth', 2, ...
            'Color', sciColors(f, :), ...
            'Marker', '.', 'MarkerSize', 15, ... % 实心点稍大一点
            'DisplayName', folderNames{f});
    end
    
    % 子图美化
    title(pointNames{p}, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 设置坐标轴属性 (加粗，向外刻度更像 SCI)
    ax = gca;
    ax.FontSize = 11;
    ax.FontName = 'Arial';
    ax.LineWidth = 1.2;     % 坐标轴线变粗
    ax.TickDir = 'out';     % 【可选】刻度朝外，看起来更干净
    
    % 统一 X 轴
    xlim([0, 2600]); 
    
    % 【修改点3】: 强制统一 Y 轴
    ylim(unified_ylim); 
    
    % 只在最左边的图显示 Y 轴标签
    if p == 1
        ylabel('Log_{10}(SNR)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        % 其他子图隐藏 Y 轴数值，显得更紧凑 (可选)
        % set(gca, 'YTickLabel', []); 
    end
    xlabel('# of basis vectors');
end

% 共享图例
lgd = legend('Orientation', 'horizontal'); 
lgd.Layout.Tile = 'south';