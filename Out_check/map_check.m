% 设置路径和序列
folder = '/data2/jiaxinli/pan/20240801_UISNR_output/20250703diag_coil/SNR_plain';  % <-- 替换为你的路径
seq = [100,200,400,800,1200,1500,2000];
slice_idx = 87;  % 你要显示的切片层

clim_linear = [0 2e-3];    % 原图色阶范围
clim_log = [-18 -11];      % log图色阶范围

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
    img_slice = squeeze(abs(img(:,:,slice_idx)));


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
        warning('文件 %s 不存在，跳过。', filepath);
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
