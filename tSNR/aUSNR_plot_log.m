function aUSNR_plot(B0,folderpath, flag, panmax, key_word, exclude_keyword)
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
    disp(['🔍 folder_tag = ' folder_tag]);
    title_prefix = ['USNR_' folder_tag ];

    fprintf('Reading data files...\n');
    for i = 1:length(panmatrix)
        key = sprintf('%s_%d', flag, panmatrix(i));
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
                data_list = zeros([data_size, 0]);
            end
            if ~isequal(size(data), data_size)
                error('Inconsistent data size: %s', file.name);
            end
            data_list(:,:,:,end+1) = data;
        end
    end

    N = size(data_list, 4);
    if N == 0
        error('No valid data files found.');
    end
    fprintf('%d data files loaded. Time elapsed: %.2f seconds\n', N, toc);

    fprintf('Building mask...\n');
    mask = all(data_list ~= 0 & ~isnan(data_list), 4);
    fprintf('Mask built. Time elapsed: %.2f seconds\n', toc);

    fprintf('Applying mask...\n');
    data_list(repmat(~mask, [1 1 1 N])) = NaN;
    fprintf('Mask applied. Time elapsed: %.2f seconds\n', toc);

    fprintf('Computing statistical maps...\n');
    mean_map = mean(data_list, 4, 'omitnan');
    std_map  = std(data_list, 0, 4, 'omitnan');
    cv_map = std_map ./ abs(mean_map);
    fprintf('Statistics computed. Time elapsed: %.2f seconds\n', toc);

    fprintf('Saving results...\n');
    save(fullfile(folderpath, sprintf('%s_mean.mat', flag)), 'mean_map');
    save(fullfile(folderpath, sprintf('%s_std.mat', flag)),  'std_map');
    save(fullfile(folderpath, sprintf('%s_cv.mat', flag)),   'cv_map');
    save(fullfile(folderpath, sprintf('%s_mask%d.mat', flag)), 'mask');
    fprintf('Results saved. Time elapsed: %.2f seconds\n', toc);

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
    vmax = -inf;
    for i = 1:N
        img = abs(get_slice(data_list(:,:,:,i)));
        data_valid = img(~isnan(img) & img ~= 0);
        if ~isempty(data_valid)
            vmax = max(vmax, max(data_valid(:)));
        end
    end

    % 👉 每个子图单独保存
    fprintf('Visualizing raw data slices (log scale)...\n');
    for i = 1:N
        fig = figure('Visible','off');
        img = log(abs(get_slice(data_list(:,:,:,i))));
        imagesc(img); axis image off;
        colormap('jet'); colorbar;
        caxis([-13 -1]);   % log形式的固定范围
        title([file_names{i} ' (log)'], 'Interpreter', 'none');
        exportgraphics(fig, fullfile(folderpath, sprintf('%s_raw_slice_log_%d.png', flag, i)), 'Resolution', 300);
        close(fig);
    end
    fprintf('Raw data visualization (log scale, separate files) completed.\n');

    fprintf('Visualizing statistical maps (log scale)...\n');
    % 平均图
    fig = figure('Visible','off'); colormap('jet');
    imagesc(log(abs(get_slice(mean_map)))); axis image off; title('Mean (log)'); colorbar;
    caxis([-13 -1]);
    exportgraphics(fig, fullfile(folderpath, sprintf('%s_mean_log.png', flag)), 'Resolution', 300);
    close(fig);

    % 标准差图
    fig = figure('Visible','off'); colormap('jet');
    imagesc(log(abs(get_slice(std_map)))); axis image off; title('Std (log)'); colorbar;
    caxis([-13 -1]);
    exportgraphics(fig, fullfile(folderpath, sprintf('%s_std_log.png', flag)), 'Resolution', 300);
    close(fig);

    % USNR图
    fig = figure('Visible','off'); colormap('jet');
    imagesc(log(abs(get_slice(cv_map)))); axis image off; title('USNR (log)'); colorbar;
    caxis([-13 -1]);
    exportgraphics(fig, fullfile(folderpath, sprintf('%s_usnr_log.png', flag)), 'Resolution', 300);
    close(fig);

    % Mask图（mask 没必要取 log，保持原样）
    fig = figure('Visible','off'); colormap('jet');
    imagesc(get_slice(mask)); axis image off; title('Mask'); colorbar;
    exportgraphics(fig, fullfile(folderpath, sprintf('%s_mask.png', flag)), 'Resolution', 300);
    close(fig);

    fprintf('Statistical visualization (log scale) completed. Total time elapsed: %.2f seconds\n', toc);
end
