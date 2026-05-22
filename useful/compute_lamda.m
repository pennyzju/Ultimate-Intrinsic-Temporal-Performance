function compute_lamda(B0,subdir,flag, seq_list)
% COMPUTE_LAMDA
% 批量处理不同 seq 的版本，计算 B1 correction 文件的变异系数 (Lamda)
%
% INPUTS:
%    B1_dir   : 包含 B1_*.mat 的目录
%    flag     : 文件夹标识 (如 'left', 'rot' 等)
%    seq_list : 序列号列表
FolderPath = B0_folderpath(B0);
B1_dir =  fullfile(FolderPath,subdir,'B1_map');
seq_list = seq_list(:);
num_seqs = numel(seq_list);

fprintf('Starting batch processing for %d sequences...\n', num_seqs);

for s = 1:num_seqs
    current_seq = seq_list(s);
    fprintf('\n[Processing seq: %d (%d/%d)]\n', current_seq, s, num_seqs);
    
    %% 1. 检查基准文件是否存在 (如果后续完全不用，此段可删除，目前仅做校验)
    baseline_name = sprintf('B1_%s_0_%d.mat', flag, current_seq);
    baseline_path = fullfile(B1_dir, baseline_name);
    
    if exist(baseline_path, 'file') ~= 2
        warning('Baseline file not found: %s.', baseline_name);
        % 如果不需要 base 叠加，这里不 continue 也可以，取决于您的逻辑需求
    end
    
    %% 2. 搜索 Correction 文件
    % 匹配模式: B1_left_*_100_correction.mat
    pattern = sprintf('B1_%s_*_%d_correction.mat', flag, current_seq);
    files = dir(fullfile(B1_dir, pattern));
    
    if isempty(files)
        warning('No correction files found for seq %d in %s.', current_seq, B1_dir);
        continue;
    end
    
    N = numel(files);
    angles = zeros(N, 1); 
    b1_all = [];
    
    %% 3. 加载 Correction 数据 (不进行叠加)
    for i = 1:N
        fname = files(i).name;
        
        % 解析文件名中的角度：匹配 B1_flag_角度_
        tok = regexp(fname, sprintf('B1_%s_(-?\\d+)_', flag), 'tokens');
        if isempty(tok)
            continue; 
        end
        angles(i) = str2double(tok{1}{1});
        
        % 加载数据
        S = load(fullfile(B1_dir, fname));
        
        % 假设 correction 文件内的变量名为 b1_map_correction
        % 如果实际变量名不同，请在此处修改
        if isfield(S, 'b1_map_correction')
            b1_data = S.b1_map_correction;
        elseif isfield(S, 'loss_map_correction') % 兼容旧名称
            b1_data = S.loss_map_correction;
        else
            % 如果不确定变量名，可以尝试取结构体第一个字段
            fields = fieldnames(S);
            b1_data = S.(fields{1});
        end
        
        % 预分配 4D 矩阵
        if isempty(b1_all)
            sz = size(b1_data);
            b1_all = zeros([sz, N], 'like', b1_data);
        end
        
        b1_all(:,:,:,i) = b1_data; % 直接赋值，不再叠加 loss0
    end
    
    %% 4. 计算 均值, 标准差 和 Lamda (CV)
    if isempty(b1_all)
        fprintf('  No valid data processed for seq %d.\n', current_seq);
        continue;
    end

    b1_mean = mean(b1_all, 4);
    b1_std  = std(b1_all, 0, 4);
    
    % 计算变异系数 Lamda
    lamda = b1_std ./ b1_mean; 
    
    fprintf('  Computed statistics over %d angles.\n', numel(angles));
    
    %% 5. 保存结果
    out_dir = B1_dir; 
    outfile = sprintf('lamda_%s_%d.mat', flag, current_seq);
    outpath = fullfile(out_dir, outfile);
    
    % 保存变量名设为 lamda
    save(outpath, 'lamda', 'b1_mean', 'b1_std', 'angles', '-v7.3');
    
    fprintf('  Saved result to: %s\n', outfile);
end

fprintf('\nBatch processing complete.\n');

end