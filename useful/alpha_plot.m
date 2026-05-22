%% 可视化多个 alpha 结果
clear
B0 = 1.5;
loss_dir = fullfile(B0_folderpath_new(B0),'\20251022seq\loss_plain'); % 替换为你的文件所在目录
flag = 'left';
seq_list = [100, 200, 400, 800, 1000, 1200, 1500, 2000,2500]; % 根据图片提供的序列
slice_idx = 48;

num_files = numel(seq_list);
rows = ceil(sqrt(num_files)); 
cols = ceil(num_files / rows);

% 预加载数据以确定全局 Colorbar 范围
all_data = cell(num_files, 1);
c_min = inf;
c_max = -inf;

fprintf('Loading data and calculating range...\n');
for i = 1:num_files
    fname = sprintf('alpha_wo_%s_%d.mat', flag, seq_list(i));
    path = fullfile(loss_dir, fname);
    
    if exist(path, 'file')
        vars = load(path);
        % 提取 (48, :, :) 切面并压缩多余维度
        slice_data = log(squeeze(vars.alpha(slice_idx, :, :)));
        all_data{i} = slice_data;
        
        % 更新全局最大最小值 (排除 NaN 和 Inf)
        valid_vals = slice_data(isfinite(slice_data));
        if ~isempty(valid_vals)
            c_min = min(c_min, min(valid_vals(:)));
            c_max = max(c_max, max(valid_vals(:)));
        end
    else
        warning('File not found: %s', fname);
    end
end

%% 开始绘制
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);
t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:num_files
    if isempty(all_data{i}), continue; end
    
    nexttile;
    % 使用 imagesc 绘制
    im = imagesc(all_data{i});
    
    % 设置统一的颜色刻度
    %clim([c_min, c_max]); 
    clim([-12, 2]); 
    colormap('jet'); % 或者使用 'hot', 'parula'
    axis image off; % 保持比例并隐藏坐标轴
    title(sprintf('Seq: %d', seq_list(i)));
end

% 在右侧添加公共 colorbar
cb = colorbar;
cb.Layout.Tile = 'east'; 
cb.Label.String = 'Alpha Value (CV)';

title(t, sprintf('Alpha Maps for Flag: %s (%dT)', flag,B0), 'FontSize', 14, 'FontWeight', 'bold');