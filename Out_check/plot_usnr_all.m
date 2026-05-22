function plot_usnr_all(folder, seq, slice_idx, line_idx, is_row)

    % 参数说明：
    % folder: 数据所在路径
    % seq: 不同序列列表，如 [100 200 400 ...]
    % slice_idx: 提取的切片层号
    % line_idx: 要提取的行/列索引
    % is_row: 是否沿行截面（true = 提取某一行，false = 提取某一列）
    
    clim_linear = [0 8e-4];    % 原图色阶范围
    clim_log = [-12 -5];      % log图色阶范围
    line_data_all = {};        % 保存所有序列的线数据
    
    % -------- 原图拼图 --------
    figure('Name','USNR Abs','Units','normalized','Position',[0.1 0.4 0.8 0.4]);
    t1 = tiledlayout(1, length(seq), 'TileSpacing','compact', 'Padding','compact');
    
    for i = 1:length(seq)
        filename = sprintf('USNR_left_0_%d.mat', seq(i));
        filepath = fullfile(folder, filename);
        if ~exist(filepath, 'file')
            warning('文件 %s 不存在，跳过。', filepath);
            continue;
        end
        data = load(filepath);
        vars = fieldnames(data);
        img = data.(vars{1});
        img_slice = abs(squeeze(img(:,:,slice_idx)));
    
        % 保存截面线数据
        if is_row
            line_data = img_slice(line_idx, 70);
        else
            line_data = img_slice(:, line_idx);
        end
        line_data_all{end+1} = line_data;
    
        % 原图显示
        nexttile(t1);
        imagesc(img_slice);
        axis image off;
        title(sprintf('seq = %d', seq(i)));
        caxis(clim_linear);
    end
    
    cb1 = colorbar;
    cb1.Location = 'eastoutside';
    cb1.Label.String = 'USNR magnitude (a.u.)';
    saveas(gcf, fullfile(folder, 'USNR_abs_all.png'));
    close;
    
    % -------- log图拼图 --------
    figure('Name','USNR Log','Units','normalized','Position',[0.1 0.4 0.8 0.4]);
    t2 = tiledlayout(1, length(seq), 'TileSpacing','compact', 'Padding','compact');
    
    for i = 1:length(seq)
        filename = sprintf('USNR_left_0_%d.mat', seq(i));
        filepath = fullfile(folder, filename);
        if ~exist(filepath, 'file')
            continue;
        end
        data = load(filepath);
        vars = fieldnames(data);
        img = data.(vars{1});
        img_slice = log(abs(squeeze(img(:,:,slice_idx))));
    
        nexttile(t2);
        imagesc(img_slice);
        axis image off;
        title(sprintf('seq = %d', seq(i)));
        caxis(clim_log);
    end
    
    cb2 = colorbar;
    cb2.Location = 'eastoutside';
    cb2.Label.String = 'log_{10}(USNR magnitude)';
    saveas(gcf, fullfile(folder, 'USNR_log_all.png'));
    close;
    
    % -------- 曲线图绘制 --------
    figure('Name','USNR Line Profile', 'Units','normalized','Position',[0.2 0.3 0.6 0.5]);
    hold on;
    for i = 1:length(line_data_all)
        plot(line_data_all{i}, 'LineWidth', 2, 'DisplayName', sprintf('seq = %d', seq(i)));
    end
    legend show;
    grid on;
    xlabel('Pixel Index');
    ylabel('USNR magnitude (a.u.)');
    title(sprintf('Line Profile at %s %d (Slice %d)', ...
        ternary(is_row,'Row','Col'), line_idx, slice_idx));
    saveas(gcf, fullfile(folder, 'USNR_line_profile.png'));
    close;
    
    end
    
    % ---- 工具函数：类似三元运算 ----
    function out = ternary(cond, a, b)
        if cond
            out = a;
        else
            out = b;
        end
    end
    