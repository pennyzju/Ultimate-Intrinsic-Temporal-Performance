function compute_alpha(B0, subdir,flag, seq_list)
% COMPUTE_ALPHA
% 批量处理不同 seq 的版本，仅计算 correction 文件中的角度
FolderPath = B0_folderpath(B0);
loss_dir =  fullfile(FolderPath,subdir,'loss_plain');

seq_list = seq_list(:);
num_seqs = numel(seq_list);

fprintf('Starting batch processing for %d sequences...\n', num_seqs);

for s = 1:num_seqs
    current_seq = seq_list(s);
    fprintf('\n[Processing seq: %d (%d/%d)]\n', current_seq, s, num_seqs);
    
    %% 1. Load baseline loss (angle = 0)
    baseline_name = sprintf('loss_%s_0_%d.mat', flag, current_seq);
    baseline_path = fullfile(loss_dir, baseline_name);
    
    if exist(baseline_path, 'file') ~= 2
        warning('Baseline file not found: %s. Skipping.', baseline_name);
        continue; 
    end
    
    S0 = load(baseline_path);
    loss0 = S0.loss_map;
    
    %% 2. Load corrected loss files
    pattern = sprintf('loss_%s_*_%d_correction.mat', flag, current_seq);
    files = dir(fullfile(loss_dir, pattern));
    
    if isempty(files)
        warning('No correction files found for seq %d.', current_seq);
        continue;
    end
    
    N = numel(files);
    % 【修改】预留 N 个位置，不再是 N+1
    angles = zeros(N, 1); 
    loss_all = [];
    
    %% 3. Load correction data
    for i = 1:N
        fname = files(i).name;
        
        % 解析角度
        tok = regexp(fname, sprintf('loss_%s_(-?\\d+)_', flag), 'tokens');
        if isempty(tok), continue; end
        angles(i) = str2double(tok{1}{1});
        
        % 加载并叠加基准值
        S = load(fullfile(loss_dir, fname));
        loss_corr = S.loss_map_correction;
        loss_sum = loss_corr + loss0;
        
        % 【修改】预分配 4D 矩阵大小为 N
        if isempty(loss_all)
            sz = size(loss_sum);
            loss_all = zeros([sz, N], 'like', loss_sum);
        end
        
        loss_all(:,:,:,i) = loss_sum;
    end
    
    %% 5. Mean, std and CV (Alpha)
    % 检查是否有有效数据，防止 mean 在空维度上操作
    if isempty(loss_all)
        continue;
    end

    loss_mean = mean(loss_all, 4);
    loss_std  = std(loss_all, 0, 4);
    
    % CV 计算
    loss_cv = loss_std ./ loss_mean; 
    
    fprintf('  Computed statistics over %d angles.\n', numel(angles));
    
    %% 6. Save results
    out_dir = loss_dir; 
    outfile = sprintf('alpha_wo_%s_%d.mat', flag, current_seq);
    outpath = fullfile(out_dir, outfile);
    
    alpha = loss_cv; %#ok<NASGU>
    save(outpath, 'alpha', 'loss_mean', 'loss_std', 'angles', '-v7.3');
    
    fprintf('  Saved result to: %s\n', outfile);
end

fprintf('\nBatch processing complete.\n');

end