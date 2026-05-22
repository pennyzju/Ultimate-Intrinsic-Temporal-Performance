function visualize_magnetic_field(M)
    % 计算磁场强度
    V = abs(M);  % 或 real(M), imag(M), angle(M)
    maxV = max(V(:));

    % 创建网格
    [x, y, z] = meshgrid(1:size(V,2), 1:size(V,1), 1:size(V,3));

    % 初始阈值（透明区域阈值）
    init_thresh = 5e-15;

    % 创建图形窗口和滑块
    fig = figure('Name','3D Magnetic Field Visualization','NumberTitle','off');
    movegui(fig, 'center');
    ax = axes('Parent', fig);
    colormap(jet); colorbar;

    % UI 滑块
    uicontrol('Style', 'text', 'Position', [20 10 140 20], ...
              'String', '透明阈值（×10^{-15}）');
    sld = uicontrol('Style', 'slider', ...
                    'Min', 0.1, 'Max', 20, 'Value', init_thresh*1e15, ...
                    'SliderStep', [0.01 0.1], ...
                    'Position', [160 10 200 20], ...
                    'Callback', @updatePlot);

    % 初始化绘图
    plotData(init_thresh);

    % 更新绘图函数
    function updatePlot(~, ~)
        thresh = sld.Value * 1e-15;  % 从滑块读取实际阈值
        cla(ax);  % 清除旧图像
        plotData(thresh);
    end

    % 主绘图函数
    function plotData(thresh)
        axes(ax);
        hold on;

        % 切片位置
        xslice = [1, round(size(V,2)/2), size(V,2)];
        yslice = [1, round(size(V,1)/2), size(V,1)];
        zslice = [1, round(size(V,3)/2), size(V,3)];

        % 切片图
        h = slice(ax, x, y, z, V, xslice, yslice, zslice);
        shading interp;

        % 设置透明度
        for i = 1:length(h)
            data = get(h(i), 'CData');
            alpha = 0.6 * ones(size(data));
            alpha(data < thresh) = 0;
            set(h(i), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 'flat', ...
                'AlphaData', alpha, ...
                'AlphaDataMapping', 'none');
        end

        % 添加等值面（isosurface）
        iso_level = thresh * 3;  % 可以根据需要调整倍数
        try
            p = patch(isosurface(x, y, z, V, iso_level));
            isonormals(x, y, z, V, p);
            set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        catch
            disp('Isosurface level太高，绘制失败（可能为空）');
        end

        camlight headlight; lighting gouraud;
        axis tight; axis equal;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('|B| with threshold = %.2e', thresh));
        view(3);
        hold off;
    end
end
