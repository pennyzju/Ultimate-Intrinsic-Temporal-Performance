%snr_process_opt_sta_seq_origin(7,'20251022seq','rot',[100,200,400,800,1000,1200,1500,2000,2200]);
%basis = [600,1400,1450,1600,1700,1800,2000];
%basis = [20:20:100];
basis=[2500];
%basis=[20:20:80]
flags = {'rot','left','front'}; % 1. 使用元胞数组 {} 定义多个字符串
B0 =1.5;
outdir = '20251022seq';
k = 2;

% 2. 使用 for 循环遍历 flags
for i = 1:numel(flags)
    curr_flag = flags{i}; % 提取当前的 flag 字符串
    
    fprintf('\n====== 正在处理 B0=%.1fT, Flag=%s ======\n', B0, curr_flag);
    
    % --- 3. 调用处理函数 ---
    % 注意：请确保你的函数内部逻辑支持传入这些参数
    %snr_process_opt_sta_seq_other(B0, outdir, curr_flag, basis);
    %snr_process_opt_sta_seq(B0, outdir, curr_flag, basis);
    % 如果 compute_alpha/lamda 内部需要路径，确保函数内已处理 B0 和 outdir 的拼接
    Angle2_compute_alpha_wo_strict_range(B0, outdir, curr_flag, basis);
%     compute_alpha_strict_range(B0, outdir, curr_flag, basis);
    Angle2_compute_lamda_strict_range(B0, outdir, curr_flag, basis);
end

fprintf('\n所有任务处理完成。\n');