function alpha_SNR_caculate_ucoilloss(B0, folderpath, flag, panmax, key_word, exclude_keyword)
    % 考虑趋肤效应增加coilloss，并计算新的SNR = B1 / sqrt(loss)
    tic;

    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test'));
    panmatrix = -panmax:1:panmax;
    file_names = {};
    data_list = [];

    %% ====== 绘制统计图 ======
    switch flag
        case {'X', 'front'}
            slice_idx = 56; get_slice = @(x) squeeze(x(:,slice_idx,:));
        case {'Y', 'left'}
            slice_idx = 48; get_slice = @(x) squeeze(x(slice_idx,:,:));
        case {'Z', 'rot'}
            slice_idx = 87; get_slice = @(x) squeeze(x(:,:,slice_idx));
        otherwise
            error('Flag must be one of: "X", "Y", "Z", "front", "left", or "rot".');
    end

    % loss 路径
    folderpath_loss = fullfile(B0_folderpath(B0), folderpath, 'loss_plain');
    % snr 路径
    folderpath_snr  = fullfile(B0_folderpath(B0), folderpath, 'SNR_plain');
    if ~exist(folderpath_snr, 'dir')
        mkdir(folderpath_snr);
    end

    folder_parts = strsplit(folderpath_loss, filesep);
    if length(folder_parts) >= 2
        folder_tag = folder_parts{end-1};
    else
        folder_tag = folder_parts{end}; 
    end
    disp(['🔍 folder_tag = ' folder_tag]);
    title_prefix = ['alpha_' folder_tag 'adducoilloss'];

    fprintf('Reading data files...\n');
    % bias 文件路径由 flag 自动决定
    bias_file = fullfile(B0_folderpath(7),'20250702dyn','loss_plain',sprintf('loss_%s_0_correction.mat', flag));
    % 安全加载
    if ~isfile(bias_file)
        error('偏移文件不存在: %s', bias_file);
    end
    biasStruct = load(bias_file);

    % 确认变量名
    if ~isfield(biasStruct, 'loss_map_correction')
        error('文件 %s 中不存在变量 loss_map_correction', bias_file);
    end
    bias = biasStruct.loss_map_correction; 

    for i = 1:length(panmatrix)
        key = sprintf('%s_%d', flag, panmatrix(i));
        if ischar(key_word)
            key_word = {key_word};
        end

        all_files = dir(fullfile(folderpath_loss, '*.mat'));
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
            filePath = fullfile(folderpath_loss, file.name);
            dataStruct = load(filePath);
            varName = fieldnames(dataStruct);
            data = dataStruct.(varName{1});
            if ~isequal(size(data), size(bias))
                error('统一偏移矩阵 bias 与数据大小不一致: %s', file.name);
            end
            data = data + sqrt(B0/7)*bias;
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

    fprintf('Computing loss statistical maps...\n');
    mean_map = mean(data_list, 4, 'omitnan');
    std_map  = std(data_list, 0, 4, 'omitnan');
    cv_map   = (std_map ./ abs(mean_map)).^2; % alpha
    fprintf('Loss statistics computed. Time elapsed: %.2f seconds\n', toc);

    fprintf('Saving loss results...\n');
    save(fullfile(folderpath_loss, sprintf('%s_mean_ucoilloss.mat', flag)), 'mean_map');
    save(fullfile(folderpath_loss, sprintf('%s_std_ucoilloss.mat', flag)),  'std_map');
    save(fullfile(folderpath_loss, sprintf('%s_cv_ucoilloss.mat', flag)),   'cv_map');
    save(fullfile(folderpath_loss, sprintf('%s_mask_ucoilloss.mat', flag)), 'mask');

    % Visualize statistics
    fprintf('Visualizing alpha statistical maps...\n');
    figure('Name', 'result', 'Position', [100, 100, 1000, 800]);
    colormap('jet');
    subplot(2,2,1); imagesc(abs(get_slice(mean_map))); axis image off; title('Mean'); colorbar;
    subplot(2,2,2); imagesc(abs(get_slice(std_map)));  axis image off; title('Std'); colorbar;
    subplot(2,2,3); imagesc(abs(get_slice(cv_map)));   axis image off; title('alpha'); colorbar;
    subplot(2,2,4); imagesc(get_slice(mask));          axis image off; title('Mask'); colorbar;
    sgtitle([title_prefix ', flag = ' flag ', adducoilloss , slice = ' num2str(slice_idx)], 'Interpreter', 'none');
    % 导出统计图
    exportgraphics(gcf, fullfile(folderpath_loss, sprintf('%s_alpha_stats_visualization_ucoilloss.png', flag)), 'Resolution', 300);

    %% ====== 新增：计算并保存 SNR map ======
    fprintf('Computing SNR maps (B1 / sqrt(loss))...\n');
    snr_list = [];

    for i = 1:N
        loss_map = data_list(:,:,:,i);

        % 找到对应的 B1 文件
        [~, fname, ~] = fileparts(file_names{i}); % e.g. loss_front_-1_correction
        angle_str = extractAfter(fname, [flag '_']); % 提取 "-1" 角度
        b1_file = fullfile(B0_folderpath(B0), folderpath, 'B1_map', ...
            sprintf('sta_b1_%s_%s.mat', flag, angle_str));

        if ~exist(b1_file, 'file')
            warning('B1 file not found: %s', b1_file);
            continue;
        end

        b1Struct = load(b1_file);
        varNameB1 = fieldnames(b1Struct);
        b1_map = b1Struct.(varNameB1{1});

        if ~isequal(size(loss_map), size(b1_map))
            error('Size mismatch between B1 and loss: %s vs %s', file_names{i}, b1_file);
        end

        snr_map =B0*B0* b1_map ./ sqrt(loss_map);
        snr_list(:,:,:,end+1) = snr_map; %#ok<AGROW>

        % 保存单个角度的 SNR
        save(fullfile(folderpath_snr, sprintf('%s_snr_ucoilloss.mat', fname)), 'snr_map');
    end

    fprintf('SNR maps computed and saved in %s.\n', folderpath_snr);

    %% ====== 计算 SNR 的 mean/std/cv ======
    fprintf('Computing SNR statistical maps...\n');
    mean_snr = mean(snr_list, 4, 'omitnan');
    std_snr  = std(snr_list, 0, 4, 'omitnan');
    cv_snr   = std_snr ./ abs(mean_snr);

    % 保存统计结果
    save(fullfile(folderpath_snr, sprintf('%s_mean_snr_ucoilloss.mat', flag)), 'mean_snr');
    save(fullfile(folderpath_snr, sprintf('%s_std_snr_ucoilloss.mat',  flag)), 'std_snr');
    save(fullfile(folderpath_snr, sprintf('%s_cv_snr_ucoilloss.mat',   flag)), 'cv_snr');

    figure('Name', 'SNR stats', 'Position', [100, 100, 1000, 800]);
    colormap('jet');
    subplot(2,2,1); imagesc(abs(get_slice(mean_snr))); axis image off; title('Mean SNR'); colorbar;
    subplot(2,2,2); imagesc(abs(get_slice(std_snr)));  axis image off; title('Std SNR'); colorbar;
    subplot(2,2,3); imagesc(abs(get_slice(cv_snr)));   axis image off; title('CV SNR'); colorbar;
    subplot(2,2,4); imagesc(get_slice(mask));          axis image off; title('Mask'); colorbar;
    sgtitle([title_prefix ', flag = ' flag ', slice = ' num2str(slice_idx)], 'Interpreter', 'none');

    exportgraphics(gcf, fullfile(folderpath_snr, sprintf('%s_snr_ucoilloss_visualization.png', flag)), 'Resolution', 300);
    close(gcf);

    fprintf('SNR statistical visualization completed. Total time elapsed: %.2f seconds\n', toc);
end
