function Angle2_compute_lamda_strict_range(B0, subdir, flag, seq_list)
% COMPUTE_LAMDA_STRICT_NO_BASE
% 1. 逻辑: 直接使用 Correction 文件计算统计量 (认为 Correction 即为 Total)
% 2. 范围: 强制 -5 到 5 度
% 3. 校验: 缺失文件报错
% 4. 处理: 剔除 0 和 NaN (Masking), 计算 Lamda = Std / Mean

    % === 参数设置 ===
    if nargin < 1, B0 = 3; end
    if nargin < 2, subdir = '20251022seq'; end
    if nargin < 3, flag = 'left'; end
    if nargin < 4, seq_list = [ 2500]; end

    angle_range = -2:1:2;
    FolderPath = B0_folderpath_new(B0);
    B1_dir = fullfile(FolderPath, subdir, 'B1_map');
    
    seq_list = seq_list(:);
    num_seqs = numel(seq_list);

    fprintf('Starting strict Lamda processing for %d sequences...\n', num_seqs);

    for s = 1:num_seqs
        current_seq = seq_list(s);
        fprintf('\n[Processing seq: %d (%d/%d)]\n', current_seq, s, num_seqs);
        
        b1_stack = [];
        
        %% 1. Loop through strict angle range
        for i = 1:length(angle_range)
            ang = angle_range(i);
            
            % 构造文件名
            corr_filename = sprintf('B1_%s_%d_%d_correction.mat', flag, ang, current_seq);
            corr_filepath = fullfile(B1_dir, corr_filename);
            
            % --- 异常检测: 文件缺失 ---
            if exist(corr_filepath, 'file') ~= 2
                error('Missing File: Angle %d file is missing for seq %d.\nPath: %s', ...
                      ang, current_seq, corr_filepath);
            end
            
            % 加载数据
            S_corr = load(corr_filepath);
            if isfield(S_corr, 'b1_map_correction')
                b1_data = S_corr.b1_map_correction;
            elseif isfield(S_corr, 'loss_map_correction')
                b1_data = S_corr.loss_map_correction;
            else
                fields = fieldnames(S_corr);
                b1_data = S_corr.(fields{1});
            end
            
            % 存入堆栈 (不叠加 Base)
            if isempty(b1_stack)
                sz = size(b1_data);
                b1_stack = zeros([sz, length(angle_range)], 'like', b1_data);
            end
            
            if ~isequal(size(b1_data), sz)
                 error('Dimension Mismatch at angle %d', ang);
            end
            
            b1_stack(:,:,:,i) = b1_data;
        end
        
        %% 2. Masking (仅剔除无效背景)
        % 即使不加 Base，去除背景的 0 值也是计算正确 Mean 的关键
        mask = all(b1_stack ~= 0 & ~isnan(b1_stack), 4);
        b1_stack(repmat(~mask, [1, 1, 1, length(angle_range)])) = NaN;
        
        %% 3. Compute Statistics
        b1_mean = mean(b1_stack, 4, 'omitnan');
        b1_std  = std(b1_stack, 0, 4, 'omitnan');
        
        % 计算 Lamda
        % 如果 b1_data 本身就是绝对值，mean 就会很大，std 相对较小
        % 得到的 Lamda 会是一个合理的百分比 (e.g., 0.05, 0.1)
        lamda = b1_std ./ abs(b1_mean); 
        
        % 清理
        lamda(b1_mean == 0) = NaN;
        
        %% 4. Save
        outfile = sprintf('Angle2_lamda_%s_%d.mat', flag, current_seq);
        outpath = fullfile(B1_dir, outfile);
        save(outpath, 'lamda', 'b1_mean', 'b1_std', 'mask', 'angle_range', '-v7.3');
        
        fprintf('  Saved result to: %s\n', outfile);
    end
    fprintf('\nBatch processing complete.\n');
end