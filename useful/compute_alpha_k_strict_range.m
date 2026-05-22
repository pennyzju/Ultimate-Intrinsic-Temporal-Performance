function compute_alpha_k_strict_range(B0, subdir, flag, seq_list,k)
% COMPUTE_ALPHA_STRICT_RANGE
% 1. 逻辑: Total = Base(0度) + Correction(theta)
% 2. 范围: 强制 -5 到 5 度
% 3. 校验: 缺失文件报错; 数据尺寸不一致报错
% 4. 处理: 剔除 0 和 NaN (Masking)

    % === 参数设置 ===
    if nargin < 1, B0 = 1.5; end
    if nargin < 2, subdir = '20251022seq'; end
    if nargin < 3, flag = 'left'; end
    %if nargin < 4, seq_list = [100,200,400,800,1000,1200,1500,2000,2500]; end
    if nargin < 4, seq_list = [2500]; end
    if nargin < 5, k = 2; end
    % 定义强制角度范围
    angle_range = -5:1:5; 

    FolderPath = B0_folderpath_new(B0);
    loss_dir = fullfile(FolderPath, subdir, 'loss_plain');
    
    seq_list = seq_list(:);
    num_seqs = numel(seq_list);

    fprintf('Starting strict batch processing for %d sequences...\n', num_seqs);
    fprintf('Angle range: %d to %d\n', min(angle_range), max(angle_range));

    for s = 1:num_seqs
        current_seq = seq_list(s);
        fprintf('\n[Processing seq: %d (%d/%d)]\n', current_seq, s, num_seqs);
        
        %% 1. Load Baseline Loss (angle = 0, plain loss)
        % 假设基准文件名为 sta_loss_flag_0.mat
        baseline_name = sprintf('loss_%s_0_%d.mat', flag,current_seq);
        baseline_path = fullfile(loss_dir, baseline_name);
        
        if exist(baseline_path, 'file') ~= 2
            error('Critical Error: Baseline file not found: %s', baseline_path);
        end
        
        S0 = load(baseline_path);
        % 这里假设变量名可能不固定，自动获取第一个变量
        varNames0 = fieldnames(S0);
        loss0 = S0.(varNames0{1}); 
        
        % 初始化数据容器
        data_stack = [];
        
        %% 2. Loop through strict angle range (-5 to 5)
        for i = 1:length(angle_range)
            ang = angle_range(i);
            
            % 构造 Correction 文件名
            % 格式假设为: sta_loss_left_-5_correction.mat
            corr_filename = sprintf('loss_%s_%d_%d_correction.mat', flag, ang,current_seq);
            corr_filepath = fullfile(loss_dir, corr_filename);
            
            % --- 异常检测 1: 文件缺失 ---
            if exist(corr_filepath, 'file') ~= 2
                error('Missing File Error: Angle %d file is missing for seq %d.\nMissing file: %s', ...
                      ang, current_seq, corr_filename);
            end
            
            % 加载 Correction 数据
            S_corr = load(corr_filepath);
            varNamesC = fieldnames(S_corr);
            loss_corr = S_corr.(varNamesC{1});
            
            % --- 异常检测 2: 尺寸不一致 ---
            if ~isequal(size(loss0), size(loss_corr))
                error('Dimension Mismatch: Angle %d correction size %s does not match baseline size %s.', ...
                      ang, mat2str(size(loss_corr)), mat2str(size(loss0)));
            end
            
            % --- 核心计算逻辑: Base + Correction ---
            total_loss = k*loss0 + loss_corr;
            
            % 存入 4D 矩阵
            if isempty(data_stack)
                % 初始化 4D 数组: [x, y, z, num_angles]
                sz = size(loss0);
                data_stack = zeros([sz, length(angle_range)], 'like', loss0);
            end
            data_stack(:,:,:,i) = total_loss;
        end
        
        %% 3. Masking & Outlier Handling (新增功能)
        fprintf('  Data loaded. Applying statistical masking...\n');
        
        % --- 异常处理 3: 构建 Mask ---
        % 规则: 只有当某像素在所有角度下都不是 0 且不是 NaN 时，才视为有效
        % dim=4 表示沿着角度维度检查
        mask = all(data_stack ~= 0 & ~isnan(data_stack), 4);
        
        % 将无效区域置为 NaN，避免污染统计结果 (比如边缘的 0 值会拉低 mean)
        data_stack(repmat(~mask, [1, 1, 1, length(angle_range)])) = NaN;
        
        %% 4. Compute Statistics
        % 使用 'omitnan' 确保即使 mask 漏掉了某些 NaN，统计函数也能正常工作
        loss_mean = mean(data_stack, 4, 'omitnan');
        loss_std  = std(data_stack, 0, 4, 'omitnan');
        
        % 计算 CV (Alpha)
        % 防止除以 0 (虽然 mask 处理过，但为了代码健壮性)
        loss_cv = (loss_std ./ abs(loss_mean)).^2;
        loss_cv(loss_mean == 0) = NaN; 
        
        fprintf('  Computed statistics over %d angles.\n', length(angle_range));
        
        %% 5. Save Results
        outfile = sprintf('alpha_%d_%s_%d.mat', k,flag, current_seq);
        outpath = fullfile(loss_dir, outfile);
        
        alpha = loss_cv; 
        % 保存 mask 和 angle_range 以便后续检查
        save(outpath, 'alpha', 'loss_mean', 'loss_std', 'mask', 'angle_range', '-v7.3');
        
        fprintf('  Saved result to: %s\n', outfile);
    end
    
    fprintf('\nStrict batch processing complete.\n');
end