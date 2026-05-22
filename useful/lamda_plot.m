%% 优化后的 Lamda 可视化 (增加基准 Mask 处理)

B0 = 1.5;
B1_dir = fullfile(B0_folderpath_new(B0),'\20251022seq\B1_map'); % 替换为你的文件所在目录
flag = 'left';
seq_list = [100, 200, 400, 800, 1000, 1200, 1500, 2000]; 
slice_idx = 48;

all_data = cell(numel(seq_list), 1);
combined_samples = []; 

%% 1. 数据加载、Mask 生成与鲁棒统计
for i = 1:numel(seq_list)
    curr_seq = seq_list(i);
    
    % --- Step A: 加载基准图 (Angle 0) 生成 Mask ---
    base_name = sprintf('B1_%s_0_%d.mat', flag, curr_seq);
    base_path = fullfile(B1_dir, base_name);
    
    if exist(base_path, 'file')
        S0 = load(base_path);
        % 假设基准变量名为 b1_map，根据实际情况修改
        base_slice = squeeze(S0.b1_map(slice_idx, :, :));
        % 生成二值 Mask：有数值的地方为 1，无数值为 0
        mask = abs(base_slice) > 0; 
    else
        warning('未找到基准文件 %s，将不使用 Mask。', base_name);
        mask = 1; % 如果没找到基准，则不遮罩
    end

    % --- Step B: 加载 Lamda 数据 ---
    fname = sprintf('lamda_%s_%d.mat', flag, curr_seq);
    path = fullfile(B1_dir, fname);
    
    if exist(path, 'file')
        vars = load(path);
        
        if isfield(vars, 'lamda')
            data_matrix = vars.lamda;
        elseif isfield(vars, 'alpha')
            data_matrix = vars.alpha;
        else
            continue;
        end
        
        % 提取切片
        slice = abs(squeeze(data_matrix(slice_idx, :, :)));
        
        % --- Step C: 应用 Mask ---
        % 将 Mask 以外的数据设为 NaN (在 imagesc 中默认显示为背景色)
        slice(~mask) = NaN; 
        all_data{i} = slice;
        
        % 收集有效数据点（仅限 Mask 内部且非 NaN 的点）
        valid_points = slice(isfinite(slice));
        combined_samples = [combined_samples; valid_points(:)];
    else
        warning('未找到文件: %s', fname);
    end
end

% 计算鲁棒的显示范围
if isempty(combined_samples)
    error('未发现有效数据点。');
end
c_min = 0; 
c_max = prctile(combined_samples, 98); 

%% 2. 绘图
figure('Color', 'w', 'Position', [100, 100, 1300, 800]);
rows = ceil(sqrt(numel(seq_list)));
cols = ceil(numel(seq_list) / rows);
t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(seq_list)
    if isempty(all_data{i}), continue; end
    
    ax = nexttile;
    % 使用 imagesc 绘制，NaN 部分会自动处理（通常显示为深蓝色或通过 'AlphaData' 隐藏）
    im = imagesc(all_data{i});
    
    % 设置背景色：让 NaN 区域显示为黑色
    set(ax, 'Color', [0 0 0]); 
    
    clim([c_min, c_max]); 
    colormap(ax, 'jet');
    axis image off;
    title(sprintf('Seq: %d', seq_list(i)));
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Lamda Value (\lambda)';

title(t, sprintf('Lamda Maps for Flag: %s (%dT)', flag,B0), 'FontSize', 14, 'FontWeight', 'bold');