% 批量读取指定文件夹中符合条件的 .mat 文件，提取某一方向上的二维切片图像，统一色阶后拼成一个总览图，并保存为图片文件。
% 输入参数：
% folderpath：字符串，指定 .mat 文件所在的文件夹路径。
% keyword：字符串，用于筛选文件名中包含该关键词的 .mat 文件。
% flag：字符串，'X' 或 'Y'，表示从三维数据中提取 X 或 Y 方向的某一切片。
% exclude_keyword（可选）：字符串或字符串 cell，表示要排除的文件名关键词。
function visualize_mat_folder(folderpath, keyword, flag, exclude_keyword)
    % 获取所有匹配的 .mat 文件
    all_files = dir(fullfile(folderpath, ['*' keyword '*.mat']));
    
    % 如果 exclude_keyword 非空，则排除包含该关键词的文件
    if nargin < 4 || isempty(exclude_keyword)
        files = all_files;
    else
        if ischar(exclude_keyword)
            exclude_keyword = {exclude_keyword};  % 转为 cell 数组
        end
        files = all_files(~contains({all_files.name}, exclude_keyword));
    end
    
    nFiles = length(files);
    if nFiles == 0
        disp('未找到符合条件的文件');
        return;
    end

    % 子图排布
    nRows = ceil(sqrt(nFiles));
    nCols = ceil(nFiles / nRows);

    % 创建 figure
    fig = figure;
    colormap('jet');
    
    % 用于 colorbar 范围统一
    clim_all = [];

    data_all = cell(nFiles, 1);
    for i = 1:nFiles
        dataStruct = load(fullfile(folderpath, files(i).name));
        varName = fieldnames(dataStruct);
        data3d = dataStruct.(varName{1});

        % 获取切片
        switch flag
            case  {'front', 'X'}
                if size(data3d, 2) < 56
                    error('数据第 %d 文件中 Y 维度小于 56', i);
                end
                slice = squeeze(data3d(:, 56, :));
            case 'Y'
                if size(data3d, 1) < 47
                    error('数据第 %d 文件中 X 维度小于 47', i);
                end
                slice = squeeze(data3d(47, :, :));
            otherwise
                error('flag 只能为 "X" 或 "Y"');
        end

        data_all{i} = slice;
        clim_all = [clim_all; abs(slice(:))];  % 取绝对值
    end

    clim = [min(clim_all), max(clim_all)];
    
    % 绘图
    for i = 1:nFiles
        subplot(nRows, nCols, i);
        imagesc(abs(data_all{i}), clim);
        axis image off;
        title(erase(files(i).name, '.mat'), 'Interpreter', 'none');
    end

    % 添加共享 colorbar
    h = colorbar('southoutside');
    h.Position = [0.35, 0.05, 0.3, 0.02];
    sgtitle(sprintf('"%s" of flag = %s ', keyword, flag));

    % 保存图片
    output_filename = fullfile(folderpath, sprintf('overview_%s_%s.png', keyword, flag));
    saveas(fig, output_filename);
    fprintf('图像已保存至：%s\n', output_filename);
end

