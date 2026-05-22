function visualize_3D_slices_title(A, mainTitle, savePath)
    % 获取维度
    [dim1, dim2, dim3] = size(A);
    x_mid = round(dim1 / 2);
    y_mid = round(dim2 / 2);
    z_mid = round(dim3 / 2);

    % 颜色范围
    vmin = min(A(:));
    vmax = max(A(:));
    threshold = 0;  % 透明阈值

    fig = figure;
    colormap(jet);

    % --- YZ 平面 ---
    ax1 = subplot(1, 3, 1);
    slice1 = permute(A(x_mid, :, :), [3, 2, 1]);
    h1 = imagesc(slice1, [vmin vmax]);
    set(h1, 'AlphaData', double(isfinite(slice1) & slice1 > threshold));
    set(gca, 'YDir', 'normal');
    hold on;
    plot([1, dim2], [z_mid, z_mid], 'w--', 'LineWidth', 1.5);
    plot([y_mid, y_mid], [1, dim3], 'w--', 'LineWidth', 1.5);
    hold off;
    axis equal off tight;
    title(sprintf('X = %d (YZ Plane)', x_mid));

    % --- XZ 平面 ---
    ax2 = subplot(1, 3, 2);
    slice2 = permute(A(:, y_mid, :), [3, 1, 2]);
    h2 = imagesc(slice2, [vmin vmax]);
    set(h2, 'AlphaData', double(isfinite(slice2) & slice2 > threshold));
    set(gca, 'YDir', 'normal');
    hold on;
    plot([1, dim1], [z_mid, z_mid], 'w--', 'LineWidth', 1.5);
    plot([x_mid, x_mid], [1, dim3], 'w--', 'LineWidth', 1.5);
    hold off;
    axis equal off tight;
    title(sprintf('Y = %d (XZ Plane)', y_mid));

    % --- XY 平面 ---
    ax3 = subplot(1, 3, 3);
    slice3 = squeeze(A(:, :, z_mid))';
    h3 = imagesc(slice3, [vmin vmax]);
    set(h3, 'AlphaData', double(isfinite(slice3) & slice3 > threshold));
    hold on;
    plot([1, dim1], [y_mid, y_mid], 'w--', 'LineWidth', 1.5);
    plot([x_mid, x_mid], [1, dim2], 'w--', 'LineWidth', 1.5);
    hold off;
    axis equal off tight;
    title(sprintf('Z = %d (XY Plane)', z_mid));

    % --- 总标题 + colorbar ---
    sgtitle(mainTitle, 'FontWeight', 'bold', 'FontSize', 14);
    cb = colorbar(ax3, 'Location', 'eastoutside');
    %cb.Label.String = '|B|';

    % --- 子图对齐 ---
    axList = [ax1, ax2, ax3];
    pos = cell2mat(get(axList, 'Position'));
    minY = min(pos(:,2));
    maxH = max(pos(:,4));
    for i = 1:3
        pos(i,2) = minY;
        pos(i,4) = maxH;
        set(axList(i), 'Position', pos(i,:));
    end

    % --- 保存为 PNG 图像（如果提供保存路径） ---
    if nargin == 3 && ~isempty(savePath)
        pngPath = fullfile(savePath, sprintf('B1_%s.png', mainTitle));
        exportgraphics(fig, pngPath, 'Resolution', 300);

        fprintf('图像已保存为：%s\n', pngPath);
    end

end
