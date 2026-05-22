clear; 

%% ================= 参数定义 =================
B0_list = [7, 3, 1.5];
a = -2;

k_list  = [49*10^(a+4), 9*10^(a+4), 1.5*1.5*10^(a+4)];

rootDir = {
    'H:\UISNR\20240801_UISNR_output\20251022seq', ...
    'G:\20251201_3T_UISNR_output\20251022seq', ...
    'F:\20251201_1p5T_UISNR_output\20251022seq'
};

flag_list = {'left','front','rot'};

% 每一列对应 idx = (b-1)*3 + f
% 每列是三张图在 subplot(3,9,?) 的位置： [SNR tSNR ratio]
% 注意：根据 sliceExpr 的顺序，实际绘图顺序是 1:tsnr, 2:SNR, 3:ratio
subplot_map = [ ...
   6  3  9  5  2  8  4  1  7;  % Row 1 (对应 tSNR)
  15 12 18 14 11 17 13 10 16;  % Row 2 (对应 SNR)
  24 21 27 23 20 26 22 19 25]; % Row 3 (对应 Ratio)

% 定义每行对应的标签名称 (对应 sliceExpr 的顺序)
row_names = {'tSNR', 'SNR', 'Ratio'};

% 每行统一 colorbar 范围（log 后）
%caxis_row = [ [-2 5]; [-2 5]; [-1 0] ];
%caxis_row = [ [2 9]; [6 13]; [-10 -1] ];
caxis_row =  [ [0 0.1]; [0 0.1]; [0.4 1] ];
%% ================= 创建 Figure =================
figure('Position',[100 100 2000 800]); 
t = tiledlayout(3,9,'TileSpacing','Compact','Padding','Compact');
% 这是一个总标题（可选）
title(t, 'Comparison of tSNR, SNR and Ratio across B0 and Orientations', 'FontSize', 16);

%% ================= 主循环 =================
for b = 1:3
    k = k_list(b);
    base = rootDir{b};
    
    % 获取当前的 B0 值，用于标签
    current_B0 = B0_list(b); 

    for f = 1:3
        flag = flag_list{f};

        % 读取数据
        SNR_data = load(fullfile(base,'SNR_plain', sprintf('USNR_%s_0_2500_correction.mat',flag)));
        lambda_data = load(fullfile(base,'B1_map', sprintf('lamda_%s_2500.mat',flag)));
        alpha_data = load(fullfile(base,'loss_plain', sprintf('alpha_%s_2500.mat',flag)));

        SNR = SNR_data.usnr_map_correction;
        lambda = lambda_data.lamda;
        alpha  = alpha_data.alpha;

        % ===== slice 定义 =====
        switch flag
            case 'left'
                sliceExpr = {
                    'squeeze(tsnr(48,14:100,13:104))', ...
                    'squeeze(k*SNR(48,14:100,13:104))', ...
                    'squeeze(Sratio(48,14:100,13:104))'};
            case 'front'
                sliceExpr = {
                    'squeeze(tsnr(14:80,56,13:104))', ...
                    'squeeze(k*SNR(14:80,56,13:104))', ...
                    'squeeze(Sratio(14:80,56,13:104))'};
            case 'rot'
                sliceExpr = {
                    'squeeze(tsnr(14:80,8:99,87))', ...
                    'squeeze(k*SNR(14:80,8:99,87))', ...
                    'squeeze(Sratio(14:80,8:99,87))'};
        end

        idx = (b-1)*3 + f;
        
        % 生成列标签文本 (例如: "7T - left")
        col_label_str = sprintf('%.1fT - %s', current_B0, flag);
        
        % 调用绘图函数，传入标签信息
        plot_tsnr_block(SNR, lambda, alpha, k, sliceExpr, subplot_map(:,idx)', caxis_row, col_label_str, row_names);
    end
end

%% ================= 函数 =================
function plot_tsnr_block(SNR, lambda, alpha, k, sliceExpr, sp_idx, caxis_row, col_str, row_strs)
% 绘制 tSNR / SNR / ratio 图块
% 增加输入参数: col_str (列名), row_strs (行名列表)

tsnr   = k * SNR ./ sqrt(1 + alpha + (lambda .* k .* SNR).^2);
Sratio = tsnr ./ (k * SNR);

data_list = {tsnr, k*SNR, Sratio};

for ii = 1:3
    nexttile(sp_idx(ii));
    img = eval(sliceExpr{ii});
    imagesc(rot90(img ));  % log 防止 0
    axis image off;
    colormap turbo;
    colorbar;
    caxis(caxis_row(ii,:));
    
    % === 修改部分：添加标签逻辑 ===
    
    % 1. 列标签 (Top Labels)
    % 逻辑：如果当前 tile 的索引 <= 9，说明它在第一行，需要加标题
    if sp_idx(ii) <= 9
        title(col_str, 'FontSize', 12, 'Interpreter', 'none');
    end
    
    % 2. 行标签 (Left Labels)
    % 逻辑：如果当前 tile 索引是 1, 10, 19 (即 mod(idx-1, 9) == 0)，说明它在最左列
    % 这里使用 text() 或者 ylabel() 均可。因为 axis off 了，ylabel 默认会隐藏。
    % 我们可以临时打开 YLabel 的显示，或者使用 text 放置在图的左侧。
    % 最简单的方法是使用 visible 的 ylabel，即使 axis off。
    
    if mod(sp_idx(ii)-1, 9) == 0
        hY = ylabel(row_strs{ii}, 'FontSize', 12, 'FontWeight', 'bold');
        % 确保 ylabel 可见，即使 axis image off
        set(hY, 'Visible', 'on'); 
    end
    
    % ==============================
end

end