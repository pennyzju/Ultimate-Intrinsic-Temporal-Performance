function snr_process_sos_dyn_Rcoilloss(B0,subdir, flag, angle_matrix)
    % folderpath: 输出根目录
    % flag: 'X'、'Y'、'Z' 或 'front'、'left'、'rot'
    % angle_deg: 最大角度（整数），会从 -angle_deg 到 angle_deg 处理
    %% 添加路径
        paths = {
            '/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'
            '/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'
            '/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'
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
    
    
    
    %% 角度列表（传入参数）；
    %angle_matrix(angle_matrix == 0) = [];
    
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
    
    %% 主循环处理每个角度
    for i = 1:length(angle_matrix)
        angle = angle_matrix(i);
        subFolder = sprintf('%s_%d', flag, angle);
        subFolderPath = fullfile(FolderPath, subFolder);
    
        if ~exist(subFolderPath, 'dir')
            fprintf('文件夹不存在，跳过: %s\n', subFolderPath);
            continue;
        end
    
        fprintf('处理子文件夹: %s\n', subFolderPath);
        matFileName = sprintf('%s_%d_material_maps.mat', flag, angle);

        try
            load(fullfile(subFolderPath, 'BASIS_B1m'),'B1m');
            load(fullfile(subFolderPath, 'object_def'),'idxS');
            load(fullfile(subFolderPath, 'BASIS_LOSS_N'),'LOSS');
            load(fullfile(subFolderPath, matFileName), 'rotated_mask');
        catch ME
            fprintf('加载基础文件时出错: %s\n', ME.message);
            continue;
        end
    
        usnr_mask_rot = usnr_mask .* rotated_mask;
        LOSS = LOSS(1:size(B1m,2), 1:size(B1m,2));

        %coilloss = 1.9117e-16;
        coilpath = fullfile(B0_folderpath(7), [flag '_0'],'BASIS_LOSS_N.mat');
        coilloss = load(coilpath);
        coilloss = sqrt(B0/7)*diag(diag(coilloss.LOSS)); 
        totLOSS = LOSS + coilloss;
    
        [usnr_map, b1_map, loss_map, rf_rate] = compute_sos_snr_basic(B1m, totLOSS, usnr_mask_rot, idxS);
    
        fprintf('当前系统时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    
        % 保存原始 map
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('USNR_%s_%d.mat', flag, angle)), 'usnr_map');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('loss_%s_%d.mat', flag, angle)), 'loss_map');
        save(fullfile(FolderPath,subdir, subfolders{3}, sprintf('mask_%s_%d.mat', flag, angle)), 'usnr_mask_rot');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('B1_%s_%d.mat', flag, angle)), 'b1_map');
        save(fullfile(FolderPath,subdir, subfolders{5}, sprintf('rfrate_%s_%d.mat', flag, angle)), 'rf_rate');
    
        % 保存角度校正 map
        [p1, p2] = define_rotation_axis(flag);
        usnr_map_correction= rotate_matrix(usnr_map, p1, p2, -angle);
        loss_map_correction= rotate_matrix(loss_map, p1, p2, -angle);
        b1_map_correction= rotate_matrix(b1_map, p1, p2, -angle);
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('USNR_%s_%d_correction.mat', flag, angle)), 'usnr_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('loss_%s_%d_correction.mat', flag, angle)), 'loss_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('b1_%s_%d_correction.mat', flag, angle)), 'b1_map_correction');
    
        clear b1_map loss_map B1m usnr_map usnr_mask_rot
        cd('..');
    end
    end
    
    