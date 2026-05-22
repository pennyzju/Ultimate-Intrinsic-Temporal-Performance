% 在指定文件夹中提取符合命名规则的一组 .mat 数据文件，对其进行统计分析（均值、标准差、变异系数），构建有效数据掩膜（mask），并对结果进行可视化与保存。
% folderpath：包含数据文件的文件夹路径。
% flag：用于构造文件名的一部分（如 'X' 或 'Y'），还用于确定切片方向。
% panmax：用于生成 -panmax:panmax 的整数列表，表示不同的偏转角度编号。
% key_word：文件名中必须包含的关键词。

% exclude_keyword（可选）：文件名中不得包含的关键词（支持单个或多个字符串）。
function alpha_caculate(folderpath, flag, panmax,key_word, exclude_keyword)
    panmatrix = -panmax:1:panmax;
    data_list = {};
    file_names = {};

    % 读取数据
    for i = 1:length(panmatrix)
        key = sprintf('%s_%d', flag, panmatrix(i));
        %all_files = dir(fullfile(folderpath, ['*' key ','keyword' *.mat']));
        all_files = dir(fullfile(folderpath, '*.mat'));
        all_files = all_files(contains({all_files.name}, key) & contains({all_files.name}, key_word));

        
        if nargin >= 4 && ~isempty(exclude_keyword)
            if ischar(exclude_keyword)
                exclude_keyword = {exclude_keyword};
            end
            all_files = all_files(~contains({all_files.name}, exclude_keyword));            
        end

        for f = 1:length(all_files)
            file = all_files(f);
            dataStruct = load(fullfile(folderpath, file.name));
            varName = fieldnames(dataStruct);
            data = dataStruct.(varName{1});
            data_list{end+1} = data;
            file_names{end+1} = file.name;
        end
    end

    if isempty(data_list)
        error('没有找到符合条件的数据文件。');
    end

    % 检查尺寸一致
    data_size = size(data_list{1});
    for i = 2:length(data_list)
        if ~isequal(size(data_list{i}), data_size)
            error('数据尺寸不一致。');
        end
    end

    % 构建 mask（所有矩阵都有效的位置）
    mask = true(data_size);
    for i = 1:length(data_list)
        %mask = mask & ~isnan(data_list{i});
        mask = mask & (data_list{i} ~= 0) & ~isnan(data_list{i});
    end

    % 合并为 4D
    N = length(data_list);
    nd = ndims(data_list{1});
    data_stack = cat(nd+1, data_list{:});
    for i = 1:N
        temp = data_stack(:,:,:,i);
        temp(~mask) = NaN;
        data_stack(:,:,:,i) = temp;
    end

    % 统计量
    mean_map = mean(data_stack, nd+1, 'omitnan');
    std_map  = std(data_stack, 0, nd+1, 'omitnan');
    cv_map   = std_map ./ (abs(mean_map) + 1e-12);

    % 保存结果
    save(fullfile(folderpath, [flag '_mean.mat']), 'mean_map');
    save(fullfile(folderpath, [flag '_std.mat']),  'std_map');
    save(fullfile(folderpath, [flag '_cv.mat']),   'cv_map');
    save(fullfile(folderpath, [flag '_mask.mat']), 'mask');

    %======== 显示切片层，根据 flag =========%
    switch upper(flag)
        case 'X'
            slice_idx = 56;  % Y = 56
            get_slice = @(x) squeeze(x(:,slice_idx,:));
            slice_dim = 2;
        case 'Y'
            slice_idx = 47;  % X = 47
            get_slice = @(x) squeeze(x(slice_idx,:,:));
            slice_dim = 1;
        otherwise
            error('flag 只能为 "X" 或 "Y"');
    end

    %======== 原始数据可视化 ========%
    %======== 原始数据可视化 ========%
    figure('Name', 'origin output + mask contour', 'Position', [100, 100, 1400, 800]);
    nCols = ceil(sqrt(N));
    nRows = ceil(N / nCols);

    for i = 1:N
        subplot(nRows, nCols, i);
        img = abs(get_slice(data_list{i}));
        
        % 设置colorbar范围仅根据数据本身的范围
        data_valid = img(~isnan(img) & img ~= 0);  % 仅保留有效数据
        caxis([min(data_valid(:)), 0.05]);  % 设置colorbar范围

        imagesc(img); axis image off;
        colorbar;
        hold on;
        contour(get_slice(mask), [0.5 0.5], 'w', 'LineWidth', 1.5);
        title(strrep(file_names{i}, '_', '\_'), 'Interpreter', 'none');
        
    end

    exportgraphics(gcf, fullfile(folderpath, [flag '_raw_with_mask.png']), 'Resolution', 300);
    % filename = flag + key_word + "_raw_with_mask.png";
    % exportgraphics(gcf, fullfile(folderpath, char(filename)), 'Resolution', 300);


    %======== 统计图可视化 ========%
    figure('Name', 'result', 'Position', [100, 100, 1000, 800]);
    colormap('jet');

    subplot(2,2,1);
    imagesc(abs(get_slice(mean_map))); axis image off; title('Mean'); colorbar;

    subplot(2,2,2);
    imagesc(abs(get_slice(std_map))); axis image off; title('Std'); colorbar;

    subplot(2,2,3);
    imagesc(abs(get_slice(cv_map))); axis image off; title('CV'); colorbar;

    subplot(2,2,4);
    imagesc(get_slice(mask)); axis image off; title('Mask'); colorbar;

    sgtitle(['flag = ' flag ', slice = ' slice_idx '）']);

    exportgraphics(gcf, fullfile(folderpath, [flag '_stats_visualization.png']), 'Resolution', 300);
    %filename = flag + key_word + "_stats_visualization.png";
    %exportgraphics(gcf, fullfile(folderpath, char(filename)), 'Resolution', 300);
end
