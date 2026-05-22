function aUSNR_check(B0,folderpath, flag, panmax, key_word, exclude_keyword)
    tic;
    panmatrix = -panmax:1:panmax;
    file_names = {};
    data_list = [];
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test'));
    folderpath = fullfile(B0_folderpath(B0),folderpath,'SNR_plain');
    
    % Extract second last folder name as a tag
    folder_parts = strsplit(folderpath, filesep);
    if length(folder_parts) >= 2
        folder_tag = folder_parts{end-1};
    else
        folder_tag = folder_parts{end};  % fallback
    end
    % 👉 Output folder_tag for verification
    disp(['🔍 folder_tag = ' folder_tag]);
    title_prefix = ['USNR_' folder_tag ];

    fprintf('Reading data files...\n');
    % Collect files matching the criteria
    for i = 1:length(panmatrix)
        key = sprintf('%s_%d', flag, panmatrix(i));
        % 支持 key_word 为字符串或字符串数组（所有关键词必须都匹配）
        if ischar(key_word)
            key_word = {key_word};
        end

        all_files = dir(fullfile(folderpath, '*.mat'));
        
        matched_flags = contains({all_files.name}, key);
        for kw = 1:length(key_word)
            matched_flags = matched_flags & contains({all_files.name}, key_word{kw});
        end

        all_files = all_files(matched_flags);

        if nargin >= 7 && ~isempty(exclude_keyword)
            if ischar(exclude_keyword)
                exclude_keyword = {exclude_keyword};
            end
            for ex = 1:length(exclude_keyword)
                all_files = all_files(~contains({all_files.name}, exclude_keyword{ex}));
            end
        end

        for f = 1:length(all_files)
            file = all_files(f);
            filePath = fullfile(folderpath, file.name);
            dataStruct = load(filePath);
            varName = fieldnames(dataStruct);
            data = dataStruct.(varName{1});
            file_names{end+1} = file.name;

            if isempty(data_list)
                data_size = size(data);
                data_list = zeros([data_size, 0]); % Initialize 4D array
            end

            if ~isequal(size(data), data_size)
                error('Inconsistent data size: %s', file.name);
            end
            data_list(:,:,:,end+1) = data; % Append as 4D array
        end
    end

    N = size(data_list, 4);
    if N == 0
        error('No valid data files found.');
    end
    fprintf('%d data files loaded. Time elapsed: %.2f seconds\n', N, toc);

    % Build a mask (positions valid in all volumes)
    fprintf('Building mask...\n');
    mask = all(data_list ~= 0 & ~isnan(data_list), 4);
    fprintf('Mask built. Time elapsed: %.2f seconds\n', toc);

    % Apply mask: set invalid regions to NaN
    fprintf('Applying mask...\n');
    data_list(repmat(~mask, [1 1 1 N])) = NaN;
    fprintf('Mask applied. Time elapsed: %.2f seconds\n', toc);

    % Compute statistics
    fprintf('Computing statistical maps...\n');
    mean_map = mean(data_list, 4, 'omitnan');
    std_map  = std(data_list, 0, 4, 'omitnan');
    cv_map = std_map ./ abs(mean_map);  % lamda
    fprintf('Statistics computed. Time elapsed: %.2f seconds\n', toc);

    % Save results
    fprintf('Saving results...\n');
    % 保存结果
    save(fullfile(folderpath, sprintf('usnr_%s_mean.mat', flag)), 'mean_map');
    save(fullfile(folderpath, sprintf('usnr_%s_std.mat', flag)),  'std_map');
    save(fullfile(folderpath, sprintf('usnr%s_cv.mat', flag)),   'cv_map');
    save(fullfile(folderpath, sprintf('usnr_%s_mask%d.mat', flag)), 'mask');
    fprintf('Results saved. Time elapsed: %.2f seconds\n', toc);

    %======= Slice visualization =========%
    switch flag
        case {'X', 'front'}
            slice_idx = 56;
            get_slice = @(x) squeeze(x(:,slice_idx,:));
        case {'Y', 'left'}
            slice_idx = 48;
            get_slice = @(x) squeeze(x(slice_idx,:,:));
        case {'Z', 'rot'}
            slice_idx = 87;
            get_slice = @(x) squeeze(x(:,:,slice_idx));
        otherwise
            error('Flag must be one of: "X", "Y", "Z", "front", "left", or "rot".');
    end

    fprintf('Visualizing raw data slices...\n');
    % Visualize original data with mask contour
    figure('Name', 'origin output + mask contour', 'Position', [100, 100, 1400, 800]);
    nCols = ceil(sqrt(N)); nRows = ceil(N / nCols);
    vmax = -inf;
    for i = 1:N
        img = abs(get_slice(data_list(:,:,:,i)));
        data_valid = img(~isnan(img) & img ~= 0);
        if ~isempty(data_valid)
            vmax = max(vmax, max(data_valid(:)));
        end
    end
    for i = 1:N
        subplot(nRows, nCols, i);
        img = abs(get_slice(data_list(:,:,:,i)));
        imagesc(img); axis image off;
        caxis([0, vmax]); 

        hold on;
        contour(get_slice(mask), [0.5 0.5], 'w', 'LineWidth', 1.5);
        title(file_names{i}, 'Interpreter', 'none');
    end

    cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);  % [left bottom width height]
    sgtitle([title_prefix ', flag = ' flag ', slice = ' num2str(slice_idx)], 'Interpreter', 'none');
    % 导出原始图像（带 mask 轮廓）
    exportgraphics(gcf, fullfile(folderpath, sprintf('usnr_%s_raw_with_mask.png', flag)), 'Resolution', 300);
    fprintf('Raw data visualization completed.\n');
    % Visualize statistics
    fprintf('Visualizing statistical maps...\n');
    figure('Name', 'result', 'Position', [100, 100, 1000, 800]);
    colormap('jet');
    subplot(2,2,1); imagesc(abs(get_slice(mean_map))); axis image off; title('Mean'); colorbar;
    subplot(2,2,2); imagesc(abs(get_slice(std_map)));  axis image off; title('Std'); colorbar;
    subplot(2,2,3); imagesc(abs(get_slice(cv_map)));   axis image off; title('USNR'); colorbar;
    subplot(2,2,4); imagesc(get_slice(mask));          axis image off; title('Mask'); colorbar;
    sgtitle([title_prefix ', flag = ' flag ', slice = ' num2str(slice_idx)], 'Interpreter', 'none');
    % 导出统计图
    exportgraphics(gcf, fullfile(folderpath, sprintf('usnr_%s_stats_visualization.png', flag)), 'Resolution', 300);
    fprintf('Statistical visualization completed. Total time elapsed: %.2f seconds\n', toc);
end
