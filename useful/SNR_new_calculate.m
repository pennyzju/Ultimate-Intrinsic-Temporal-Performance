%% 批量计算 SNR (B1 / sqrt(loss + loss0)) - 支持多序列循环
clear; clc; close all;

%% ================= 1. 参数设置 =================
% 循环列表
angle_list = 0;           % 角度: -5 到 5
B0_list    = [1.5, 3, 7];      % 场强
flag_list  = {'left', 'front', 'rot'}; 
seq_list   = [20:20:100,200, 400, 800, 1000, 1200, 1500, 2000, 2200,2500];           % [新增] 序列列表，例如 [1000, 2500]

% 子目录名称
subdir = '20251022seq'; 

%% ================= 2. 主处理循环 =================
for i = 1:length(B0_list)
    curr_B0 = B0_list(i);
    
    % 获取当前 B0 的根目录
    try
        base_path = B0_folderpath_new(curr_B0); 
    catch
        error('未找到 B0_folderpath_new 函数，请确保其在路径中或手动设置 base_path');
    end
    
    % 构建输入/输出路径
    B1_dir   = fullfile(base_path, subdir, 'B1_map');
    loss_dir = fullfile(base_path, subdir, 'loss_plain');
    out_dir  = fullfile(base_path, subdir, 'SNR_plain_new');
    
    % 如果输出目录不存在，则创建
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
        fprintf('创建输出目录: %s\n', out_dir);
    end
    
    fprintf('正在处理 B0 = %.1f T ...\n', curr_B0);
    
    % --- Loop 2: 序列循环 (新增) ---
    for s = 1:length(seq_list)
        curr_seq = seq_list(s);
        fprintf('  Sequence: %d\n', curr_seq);
        
        % --- Loop 3: Flag 循环 ---
        for j = 1:length(flag_list)
            curr_flag = flag_list{j};
            
            % =======================================================
            % [关键] 加载 loss0 (基准损耗)
            % 注意：loss0 依赖于 flag 和 seq，所以放在这里加载
            % =======================================================
            try
                % 假设 baseline 文件名中的 angle 固定为 0
                baseline_name = sprintf('loss_%s_0_%d.mat', curr_flag, curr_seq);
                
                % [注意] 您代码中指定从 7T 的文件夹读取 loss0，此处保留该逻辑
                % 如果 loss0 应该随 curr_B0 变化，请将下面的 B0_folderpath_new(7) 改为 base_path
                baseline_path = fullfile(B0_folderpath_new(7), subdir, 'loss_plain', baseline_name);
                
                if exist(baseline_path, 'file')
                    S0 = load(baseline_path);
                    varNames0 = fieldnames(S0);
                    loss0 = S0.(varNames0{1}); 
                else
                    warning('    [警告] loss0 文件缺失，跳过当前 Flag: %s', baseline_path);
                    continue; % 如果没有 loss0，跳过当前 flag 的所有角度
                end
            catch ME
                warning('    [错误] 读取 loss0 失败: %s', ME.message);
                continue;
            end
            % =======================================================
            
            % --- Loop 4: 角度循环 ---
            for k = 1:length(angle_list)
                curr_angle = angle_list(k);
                
                % --- 2.1 构建文件名 (将 2500 替换为 %d) ---
                % 规则: B1_flag_angle_seq_correction.mat
                file_B1   = sprintf('B1_%s_%d_%d.mat', curr_flag, curr_angle, curr_seq);
                file_loss = sprintf('loss_%s_%d_%d.mat', curr_flag, curr_angle, curr_seq);
                file_out  = sprintf('USNR_%s_%d_%d.mat', curr_flag, curr_angle, curr_seq);
                
                path_B1   = fullfile(B1_dir, file_B1);
                path_loss = fullfile(loss_dir, file_loss);
                path_out  = fullfile(out_dir, file_out);
                
                % --- 2.2 检查文件是否存在 ---
                if exist(path_B1, 'file') && exist(path_loss, 'file')
                    
                    % --- 2.3 读取数据 ---
                    data_B1_struct = load(path_B1);
                    data_loss_struct = load(path_loss);
                    
                    % 提取 B1 矩阵
                    field_names_b1 = fieldnames(data_B1_struct);
                    val_B1 = data_B1_struct.(field_names_b1{1}); 
                    
                    % 提取 Loss 矩阵
                    field_names_loss = fieldnames(data_loss_struct);
                    val_loss = data_loss_struct.(field_names_loss{1});
                    
                    % --- 2.4 计算 SNR ---
                    if ~isequal(size(val_B1), size(val_loss))
                        warning('B1 和 Loss 维度不匹配: %s', file_B1);
                        continue;
                    end
                    
                    try
                        % 公式: SNR = B1 / sqrt(loss + loss0)
                        % 确保 loss0 维度与 val_loss 兼容 (标量或同尺寸矩阵)
                        snr = abs(val_B1) ./ sqrt(abs(val_loss) + abs(loss0)); 
                    catch ME
                        warning('计算出错 (%s): %s', file_B1, ME.message);
                        continue;
                    end
                    
                    snr(~isfinite(snr)) = 0;
                    
                    % --- 2.5 保存结果 ---
                    save(path_out, 'snr');
                    
                else
                    % fprintf('    [跳过] 文件缺失: %s \n', file_B1);
                end
            end % End Angle
        end % End Flag
    end % End Seq
end % End B0

fprintf('全部处理完成！\n');