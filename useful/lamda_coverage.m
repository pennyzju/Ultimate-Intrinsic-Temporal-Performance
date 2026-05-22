%% Lamda 均值收敛曲线绘制 (小 Mask 区域)
B1_dir = 'H:\UISNR\20240801_UISNR_output\20251022seq\B1_map'; 
flag = 'left';
% 包含更多采样点以观察数值收敛过程
seq_list = [100, 200, 400, 600, 800, 1000, 1200, 1400, 1450, 1500, 1600, 1700, 1800, 2000];
slice_idx = 48; 

mean_values = zeros(size(seq_list)); 

% 定义收缩结构元素 (2像素半径)
se = strel('disk', 2); 

for i = 1:numel(seq_list)
    fname = sprintf('lamda_%s_%d.mat', flag, seq_list(i));
    path = fullfile(B1_dir, fname);
    
    if exist(path, 'file')
        data = load(path);
        
        % 提取变量
        if isfield(data, 'lamda')
            temp_lamda = squeeze(data.lamda(slice_idx, :, :));
        elseif isfield(data, 'alpha')
            temp_lamda = squeeze(data.alpha(slice_idx, :, :));
        else
            warning('文件 %s 中未找到变量', fname);
            mean_values(i) = NaN;
            continue;
        end
        
        % --- 关键：实现小 Mask 均值计算 ---
        % 1. 生成初始二值 Mask (排除 0 值)
        initial_mask = abs(temp_lamda) > 0;
        
        if any(initial_mask(:))
            % 2. 填补内部孔洞，确保腐蚀只发生在外部边缘
            solid_mask = imfill(initial_mask, 'holes');
            
            % 3. 腐蚀 Mask，得到缩小 2 像素后的区域
            shrunk_mask = imerode(solid_mask, se);
            
            % 4. 结合有效值检查 (排除 NaN/Inf)
            valid_roi = shrunk_mask & isfinite(temp_lamda);
            
            if any(valid_roi(:))
                % 5. 仅计算收缩区域内的均值 (取绝对值处理复数)
                mean_values(i) = mean(abs(temp_lamda(valid_roi)));
            else
                mean_values(i) = NaN;
            end
        else
            mean_values(i) = NaN;
        end
    else
        mean_values(i) = NaN;
    end
end

%% 绘图部分
figure('Color', 'w', 'Name', 'Lamda Convergence (Small Mask)');
plot(seq_list, mean_values, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', '#D95319', 'Color', '#D95319');

grid on;
xlabel('Number of Basis Vectors (seq)', 'FontSize', 12);
ylabel('Mean Lamda (\lambda) in Shrunk ROI', 'FontSize', 12);
title(sprintf('Numerical Convergence of \\lambda (Slice %d, 2px Eroded Mask)', slice_idx), 'FontSize', 14);

% 绘制收敛水平线
last_val = mean_values(find(~isnan(mean_values), 1, 'last'));
if ~isempty(last_val)
    yline(last_val, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Plateau');
end