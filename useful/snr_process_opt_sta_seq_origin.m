function snr_process_opt_sta_seq_origin(B0, subdir, flag, seqlist)
    % folderpath: 输出根目录
    % flag: 'X'、'Y'、'Z' 或 'front'、'left'、'rot'
    % angle_deg: 最大角度（整数），会从 -angle_deg 到 angle_deg 处理
    %对角线上增加7T的noise-cov
    %% 添加路径
        paths = {
            'E:\UISNR/MARIE_ubasis_v6_noSVD_BIGMEM'
            'E:\UISNR/unaccelerated_matlab'
            'E:\UISNR/RUN_DUKE_BASIS_SET'
        };
        
        for i = 1:length(paths)
            addpath(genpath(paths{i}));
        end
    
    %构建子文件夹
    FolderPath = B0_folderpath(B0);

    subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map','rf_map'};
    for i = 1:length(subfolders)
        folder_name = fullfile(FolderPath,subdir, subfolders{i});
        if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    end
    
    %加载所有mask合并组成的bigmask
    mask_path = fullfile(FolderPath, 'left_0', 'object_def_bigmask.mat');
    load(mask_path, 'maskall','idxS');
    
    %% 创建 mask
    mask_filename = ['origin_' flag '.mat'];
    mask_path = fullfile(FolderPath,subdir, subfolders{3}, mask_filename);
    
    if exist(mask_path, 'file')
        fprintf('文件已存在，直接加载：%s\n', mask_path);
        load(mask_path); % 加载 usnr_mask
    else
        usnr_mask = zeros(95,112,117);
        switch flag
            case {'front', 'X'}
                usnr_mask(:,56,:) = 1;
            case {'left', 'Z'}
                usnr_mask(48,:,:) = 1;
            case {'rot', 'Y'}
                usnr_mask(:,:,87) = 1;
            otherwise
                error('未知的 flag 值，请使用 "front"、"left"、"rot" 或 "X"、"Y"、"Z"');
        end
        save(mask_path, 'usnr_mask');
        fprintf('已生成并保存 usnr_mask：%s\n', mask_path);
    end

    usnr_mask_seq = usnr_mask .* maskall;

    %初始位置多位置下的基函数集
    subFolderPath = fullfile(FolderPath, 'left_0','bigmask');
    
    if ~exist(subFolderPath, 'dir')
        fprintf('文件夹不存在，跳过: %s\n', subFolderPath);
    end
        
    fprintf('处理子文件夹: %s\n', subFolderPath);
    %matFileName = sprintf('%s_%d_material_maps.mat', flag, angle);

    try
        load(fullfile(subFolderPath, 'BASIS_B1m'),'B1m');
        load(fullfile(subFolderPath, 'BASIS_LOSS_N'),'LOSS');
    catch ME
        fprintf('加载基础文件时出错: %s\n', ME.message);
    end
        LOSS = LOSS(1:size(B1m,2), 1:size(B1m,2));
        
    angle = 0;

    for seq = seqlist
        B1mseq = B1m(:,1:seq);
        LOSSseq = LOSS(1:seq,1:seq);
        [usnr_map, b1_map, loss_map, rf_rate, rf_map] = compute_usnr_unacc_fast_v5(B1mseq, LOSSseq, usnr_mask_seq, idxS);

        fprintf('当前系统时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

        % 保存原始 map
        save(fullfile(FolderPath, subdir, subfolders{1}, sprintf('USNR_%s_%d_%d.mat', flag, angle, seq)), 'usnr_map');
        save(fullfile(FolderPath, subdir, subfolders{2}, sprintf('loss_%s_%d_%d.mat', flag, angle, seq)), 'loss_map');
        save(fullfile(FolderPath, subdir, subfolders{3}, sprintf('mask_%s_%d_%d.mat', flag, angle, seq)), 'usnr_mask_seq');
        save(fullfile(FolderPath, subdir, subfolders{4}, sprintf('B1_%s_%d_%d.mat', flag, angle, seq)), 'b1_map');
        save(fullfile(FolderPath, subdir, subfolders{5}, sprintf('rf_%s_%d_%d.mat', flag, angle, seq)), 'rf_map');
        save(fullfile(FolderPath, subdir, subfolders{5}, sprintf('rfrate_%s_%d_%d.mat', flag, angle, seq)), 'rf_rate');

    end
    clear b1_map loss_map B1m usnr_map usnr_mask_rot
    cd('..');
end

    
    