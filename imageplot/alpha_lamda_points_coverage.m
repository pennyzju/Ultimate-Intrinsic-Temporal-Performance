%% 多场强 CV 收敛对比（方案 B：固定点，点级收敛）
clear;

%% ================= 基本参数 =================
% 请确保 B0_folderpath 函数在您的路径中可用
B0_list = [1.5, 3, 7];
subdir  = '20251022seq';
flag    = 'left';

% 序列列表
seq_list = [ 100, 200, 400, 800, 1000, 1200, 1500, 2000, 2200, 2500];

%% ================= 方案 B：人工指定代表性点 =================
% 每一行是一个体素坐标 [x y z]
point_list = [
    48 56 74;   % 中心
    %48 60 74;   % 右
    48 70 74; 
    %48 80 74;
    48 84 74;% 左
];

num_points = size(point_list,1);

%% ================= Step 1：生成线性索引（只做一次） =================
test_B0  = B0_list(1);
test_seq = seq_list(1);

% 确保 B0_folderpath 函数存在，这里假设它能正常工作
try
    test_dir = fullfile(B0_folderpath_new(test_B0), subdir, 'loss_plain');
catch
    error('请确保 B0_folderpath 函数在 MATLAB 路径中。');
end

test_file = fullfile(test_dir, sprintf('alpha_wo_%s_%d.mat', flag, test_seq));

assert(exist(test_file,'file')==2, 'Reference file not found: %s', test_file);

tmp = load(test_file);
if isfield(tmp,'alpha')
    ref = tmp.alpha;
else
    ref = tmp.lamda;
end

data_size = size(ref);

% 坐标合法性检查
assert(all(point_list(:,1) <= data_size(1)) && ...
       all(point_list(:,2) <= data_size(2)) && ...
       all(point_list(:,3) <= data_size(3)), ...
       'Point coordinates exceed data size.');

% 转为线性索引
sample_idx = sub2ind(data_size, ...
    point_list(:,1), point_list(:,2), point_list(:,3));

%% ================= 预分配（B0 × point × seq） =================
all_alpha = nan(numel(B0_list), num_points, numel(seq_list));
all_lamda = nan(numel(B0_list), num_points, numel(seq_list));

%% ================= Step 2：数据收集（逐点） =================
for b = 1:numel(B0_list)
    curr_B0 = B0_list(b);

    loss_dir = fullfile(B0_folderpath(curr_B0), subdir, 'loss_plain');
    B1_dir   = fullfile(B0_folderpath(curr_B0), subdir, 'B1_map');

    fprintf('Processing B0 = %.1f T...\n', curr_B0);

    for i = 1:numel(seq_list)
        curr_seq = seq_list(i);

        %% -------- Alpha (Loss CV) --------
        alpha_file = fullfile(loss_dir, ...
            sprintf('alpha_%s_%d.mat', flag, curr_seq));

        if exist(alpha_file,'file')
            a_data = load(alpha_file);
            if isfield(a_data,'alpha')
                temp_a = abs(a_data.alpha);
            else
                temp_a = abs(a_data.lamda);
            end

            vals = temp_a(sample_idx);
            vals(~isfinite(vals) | vals <= 0) = NaN;

            for p = 1:num_points
                all_alpha(b,p,i) = vals(p);
            end
        end

        %% -------- Lamda (B1 CV) --------
        lamda_file = fullfile(B1_dir, ...
            sprintf('lamda_%s_%d.mat', flag, curr_seq));

        if exist(lamda_file,'file')
            l_data = load(lamda_file);
            if isfield(l_data,'lamda')
                temp_l = abs(l_data.lamda);
            else
                temp_l = abs(l_data.alpha);
            end

            vals = temp_l(sample_idx);
            vals(~isfinite(vals) | vals <= 0) = NaN;

            for p = 1:num_points
                all_lamda(b,p,i) = vals(p);
            end
        end
    end
end

fprintf('Data processing complete.\n');

%% ================= Step 3：可视化绘图 =================

%% part 3.1: Alpha (Loss CV) 绘图
figure('Color','w','Position',[100 500 1200 400], 'Name', 'Alpha Convergence');
t1 = tiledlayout(1, numel(B0_list), 'TileSpacing','compact','Padding','compact');

% --- [关键] 计算 Alpha 全局统一 Y 轴范围 ---
valid_data_a = all_alpha(~isnan(all_alpha));
if isempty(valid_data_a)
    ylim_alpha = [0 1];
else
    y_min = min(valid_data_a);
    y_max = max(valid_data_a);
    margin = (y_max - y_min) * 0.1; 
    if margin == 0, margin = 0.1 * abs(y_max); end
    if margin == 0, margin = 0.1; end % 防止全是0的情况
    ylim_alpha = [y_min - margin, y_max + margin];
end
% -------------------------------------

for b = 1:numel(B0_list)
    ax = nexttile(t1); hold on;

    for p = 1:num_points
        % 绘制曲线 (实心圆点)
        pl = plot(seq_list, squeeze(all_alpha(b,p,:)), ...
            '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
        set(pl, 'MarkerFaceColor', pl.Color); 
    end

    % --- 样式调整 ---
    grid off;      % 取消网格
    box off;       % 取消边框
    
    xlabel('Sequence Index');
    ylabel('\alpha (Loss CV)');
    title(sprintf('B_0 = %.1f T', B0_list(b)));
    
    % 设置统一纵坐标
    ylim(ylim_alpha); 
    
    % 美化坐标轴
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11);

    % 图例仅在第一个子图显示
    if b == 1
        legend(arrayfun(@(x) sprintf('Point %d', x), ...
            1:num_points, 'UniformOutput', false), ...
            'Location','northeast', 'Box', 'off'); 
    end
end

title(t1, ...
    sprintf('\\alpha Convergence at Fixed Spatial Points (N = %d)', num_points), ...
    'FontSize',14,'FontWeight','bold');



%% part 3.2: Lamda (B1 CV) 绘图 (新增部分)
figure('Color','w','Position',[100 100 1200 400], 'Name', 'Lamda Convergence');
t2 = tiledlayout(1, numel(B0_list), 'TileSpacing','compact','Padding','compact');

% --- [关键] 计算 Lamda 全局统一 Y 轴范围 ---
valid_data_l = all_lamda(~isnan(all_lamda));
if isempty(valid_data_l)
    ylim_lamda = [0 1];
else
    y_min = min(valid_data_l);
    y_max = max(valid_data_l);
    margin = (y_max - y_min) * 0.1; 
    if margin == 0, margin = 0.1 * abs(y_max); end
    if margin == 0, margin = 0.1; end
    ylim_lamda = [y_min - margin, y_max + margin];
end
% -------------------------------------

for b = 1:numel(B0_list)
    ax = nexttile(t2); hold on;

    for p = 1:num_points
        % 绘制曲线 (实心圆点)
        pl = plot(seq_list, squeeze(all_lamda(b,p,:)), ...
            '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
        set(pl, 'MarkerFaceColor', pl.Color); % 实心点设置
    end

    % --- 样式调整 (与 Alpha 保持一致) ---
    grid off;      % 取消网格
    box off;       % 取消边框
    
    xlabel('Sequence Index');
    ylabel('\lambda (B_1 CV)'); % 更改 Y 轴标签
    title(sprintf('B_0 = %.1f T', B0_list(b)));
    
    % 设置统一纵坐标 (使用 Lamda 的范围)
    ylim(ylim_lamda); 
    
    % 美化坐标轴
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11);

    % 图例仅在第一个子图显示
    if b == 1
        legend(arrayfun(@(x) sprintf('Point %d', x), ...
            1:num_points, 'UniformOutput', false), ...
            'Location','northeast', 'Box', 'off'); 
    end
end

title(t2, ...
    sprintf('\\lambda Convergence at Fixed Spatial Points (N = %d)', num_points), ...
    'FontSize',14,'FontWeight','bold');