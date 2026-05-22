%% ================= 双K值对比绘图 (3x6 布局) =================
clear; clc;

%% 1. 基础参数定义
k_list_1 = [1.5*1.5*1e3, 9*1e3, 49*1e3];
k_list_2 = [1.5*1.5*1e6, 9*1e6, 49*1e6];

configs(1).name = '1.5 T';
configs(1).path = 'F:\20251201_1p5T_UISNR_output\20251022seq';
configs(1).idx  = 1;

configs(2).name = '3 T';
configs(2).path = 'G:\20251201_3T_UISNR_output\20251022seq';
configs(2).idx  = 2;

configs(3).name = '7 T';
configs(3).path = 'H:\UISNR\20240801_UISNR_output\20251022seq';
configs(3).idx  = 3;

flag = 'left';

% 左侧范围
clim_left.tsnr  = [-4, 1];
clim_left.snr   = [-4, 1];
clim_left.ratio = [0.4, 1];

% 右侧范围
clim_right.tsnr  = [0.5, 5];
clim_right.snr   = [0.5, 5];
clim_right.ratio = [0, 1];

%% 2. 初始化画布
figure('Color','w','Position',[50 50 1600 800]);
t = tiledlayout(3,6,'TileSpacing','compact','Padding','compact');
title(t,'Comparison of tSNR, SNR, and Ratio (Sagittal View)', ...
      'FontSize',16,'FontWeight','bold');

row_labels = {'log_{10}(tSNR)','log_{10}(SNR)','Ratio (tSNR/SNR)'};

%% 3. 主循环
for col = 1:6
    
    % ===== 判断左右组 =====
    if col <= 3
        field_idx = col;
        curr_k = k_list_1(configs(field_idx).idx);
        curr_clim = clim_left;
        group_name = 'Condition A';
    else
        field_idx = col - 3;
        curr_k = k_list_2(configs(field_idx).idx);
        curr_clim = clim_right;
        group_name = 'Condition B';
    end
    
    curr_config = configs(field_idx);
    base_path = curr_config.path;
    
    % ===== 读取数据 =====
    SNR_data    = load(fullfile(base_path,'SNR_plain', ...
                    sprintf('USNR_%s_0_2500_correction.mat',flag)));
    lambda_data = load(fullfile(base_path,'B1_map', ...
                    sprintf('lamda_%s_2500.mat',flag)));
    alpha_data  = load(fullfile(base_path,'loss_plain', ...
                    sprintf('alpha_%s_2500.mat',flag)));

    SNR   = SNR_data.usnr_map_correction;
    lamda = lambda_data.lamda;
    alpha = alpha_data.alpha;

    % ===== 将 (:,87,:) 清零 =====
    if ndims(SNR)==3,   SNR(:,:,17)   = 0; end
    if ndims(lamda)==3, lamda(:,:,17)= 0; end
    if ndims(alpha)==3, alpha(:,:,17)= 0; end

    % ===== ROI =====
    slice_SNR = get_roi_data(SNR,flag);
    slice_lam = get_roi_data(lamda,flag);
    slice_alp = get_roi_data(alpha,flag);

    % ===== 计算物理量 =====
    img_SNR = curr_k .* slice_SNR;
    denom   = sqrt(1 + slice_alp + (slice_lam .* img_SNR).^2);
    img_tSNR = img_SNR ./ denom;
    img_Ratio = img_tSNR ./ img_SNR;

    % ===== Row 1: tSNR =====
    nexttile((0*6)+col);
    plot_slice(log10(img_tSNR), curr_clim.tsnr);
    if col==1
        ylabel(row_labels{1},'FontSize',12,'FontWeight','bold');
    end
    title(sprintf('%s\n%s',group_name,curr_config.name),'FontSize',11);

    % ===== Row 2: SNR =====
    nexttile((1*6)+col);
    plot_slice(log10(img_SNR), curr_clim.snr);
    if col==1
        ylabel(row_labels{2},'FontSize',12,'FontWeight','bold');
    end

    % ===== Row 3: Ratio =====
    nexttile((2*6)+col);
    plot_slice(img_Ratio, curr_clim.ratio);
    if col==1
        ylabel(row_labels{3},'FontSize',12,'FontWeight','bold');
    end
end


%% 4. 每三个子图共享一个 colorbar（稳定版）

for row = 1:3
    
    % ==== 左组三个子图 ====
    ax_left = [];
    for c = 1:3
        ax_left = [ax_left, nexttile((row-1)*6 + c)];
    end
    cb_left = colorbar(ax_left(end));
    cb_left.Location = 'eastoutside';
    cb_left.FontSize = 10;
    
    % ==== 右组三个子图 ====
    ax_right = [];
    for c = 4:6
        ax_right = [ax_right, nexttile((row-1)*6 + c)];
    end
    cb_right = colorbar(ax_right(end));
    cb_right.Location = 'eastoutside';
    cb_right.FontSize = 10;
    
end


%% 5. 中间分割线
annotation('line',[0.5 0.5],[0.05 0.92], ...
           'Color','k','LineStyle','--','LineWidth',1.5);

%% ================= 辅助函数 =================
function plot_slice(img, clim_val)
    imagesc(rot90(img));
    axis image off;
    colormap(gca,'turbo');
    clim(clim_val);
end

function slice_out = get_roi_data(data, flag)

    if ndims(data)==3
        switch flag
            case 'left'
                slice_out = squeeze(data(48,14:100,13:104));
            case 'front'
                slice_out = squeeze(data(14:80,56,13:104));
            case 'rot'
                slice_out = squeeze(data(14:80,8:99,87));
            otherwise
                slice_out = data;
        end
    else
        slice_out = data;
    end
end