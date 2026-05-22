function rotate_snr_process_mag(FolderPath, flag, angle_matrix,subdir, abs_rf,abs_q)
    subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map'};
    
    for i = 1:length(subfolders)
        folder_name = fullfile(FolderPath,subdir, subfolders{i});
        if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    end

    add_custom_paths();

    %angle_matrix = -angle_deg:1:angle_deg;
    %angle_matrix = angle_deg;
    %创建或加载原始 mask
    usnr_mask = prepare_origin_mask(FolderPath, subdir, subfolders{3}, flag);

    for angle = angle_matrix
        subfolder = sprintf('%s_%d', flag, angle);
        subFolderPath = fullfile(FolderPath, subfolder);
        
        if ~exist(subFolderPath, 'dir')
            fprintf('跳过不存在的子文件夹: %s\n', subFolderPath);
            continue;
        end
        fprintf('处理子文件夹: %s\n', subFolderPath);

        % 查找并加载 material 文件
        matFile = sprintf('%s_%d_material_maps.mat', flag, angle);
        
        try
            load(fullfile(subFolderPath, matFile));  % 加载 pan_mask
        catch ME
            fprintf('加载文件 %s 出错: %s\n', matFile, ME.message); 
            continue;
        end

        % 加载 object_def.mat
        try
            load(fullfile(subFolderPath, 'object_def.mat'));
        catch ME
            fprintf('加载 object_def.mat 出错: %s\n', ME.message); 
            continue;
        end

        

        try
            magFolder = fullfile(subFolderPath, 'mag');
            data = load_all_files(magFolder);
        catch ME
            fprintf('加载数据失败: %s\n', ME.message);
            continue;
        end

        usnr_mask_rot = usnr_mask .* rotated_mask;
        LOSS = data.LOSS(1:size(data.B1m,2), 1:size(data.B1m,2));

        %[usnr_map, b1_map, loss_map] = compute_usnr_unacc_fast_v4(data.B1m, LOSS, usnr_mask_rot, data.idxS);
        % 传入计算函数，并传递 abs_rf 和 abs_q 参数
        [usnr_map, b1_map, loss_map, rf_rate] = compute_sos_snr_general(data.B1m, LOSS, usnr_mask_rot, idxS, abs_rf, abs_q);
        fprintf('时间戳: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));


        % 保存原始 map
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('USNR_%s_%d.mat', flag, angle)), 'usnr_map');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('loss_%s_%d.mat', flag, angle)), 'loss_map');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('B1_%s_%d.mat', flag, angle)), 'b1_map');

        % 保存角度校正 map
        [p1, p2] = define_rotation_axis(flag);
        usnr_map_correction= rotate_matrix(usnr_map, p1, p2, -angle);
        loss_map_correction= rotate_matrix(loss_map, p1, p2, -angle);
        b1_map_correction= rotate_matrix(b1_map, p1, p2, -angle);
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('USNR_%s_%d_correction.mat', flag, angle)), 'usnr_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('loss_%s_%d_correction.mat', flag, angle)), 'loss_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('b1_%s_%d_correction.mat', flag, angle)), 'b1_map_correction');
        
        clear usnr_map b1_map loss_map usnr_mask_rot
        cd('..');
    end
end

function add_custom_paths()
    paths = {
        '/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'
        '/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'
        '/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'
    };
    for i = 1:length(paths)
        addpath(genpath(paths{i}));
    end
end

function usnr_mask = prepare_origin_mask(basePath, subdir, maskFolder, flag)
    mask_path = fullfile(basePath, subdir, maskFolder, ['origin_' flag '.mat']);
    if exist(mask_path, 'file')
        load(mask_path, 'usnr_mask');
        return;
    end

    usnr_mask = zeros(95,112,117);
    switch flag
        case 'front', usnr_mask(:,56,:) = 1;
        case 'left',  usnr_mask(48,:,:) = 1;
        case 'rot',   usnr_mask(:,:,87) = 1;
        otherwise, error('flag 仅支持 front / left / rot');
    end
    save(mask_path, 'usnr_mask');
end

function data = load_all_files(FolderPath)
    cd(FolderPath);
    load('BASIS_B1m.mat', 'B1m');
    load('BASIS_LOSS_N.mat', 'LOSS');

    data.B1m = B1m;
    data.LOSS = LOSS;
end

function [p1, p2] = define_rotation_axis(flag)
    m = 95; n = 112; p = 117;
    switch flag
        case 'left'
            p1 = [0, n/2, 30]; p2 = [m, n/2, 30];
        case 'front'
            p1 = [m/2, 0, 30]; p2 = [m/2, n, 30];
        case 'rot'
            p1 = [m/2, n/2, 0]; p2 = [m/2, n/2, p];
        otherwise
            error('未知的 flag 值');
    end
end

