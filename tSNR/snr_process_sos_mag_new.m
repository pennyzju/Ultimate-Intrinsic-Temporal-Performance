function snr_process_sos_mag_new(FolderPath, flag, angle_deg, subdir,abs_rf, abs_q)

    disp(['FolderPath = ', FolderPath]);
    subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map','rf_map'};
    for i = 1:length(subfolders)
        folder_name = fullfile(FolderPath, subdir, subfolders{i});
        if exist(folder_name, 'dir')
            if is_writable(folder_name)
                disp(['Folder exists and is writable: ', folder_name]);
            else
                warning(['Folder exists but not writable: ', folder_name]);
            end
        else
            try
                mkdir(folder_name);
            catch ME
                warning('Failed to create folder: %s', ME.message);
            end
        end
    end

    %% 添加路径（你的路径根据需要修改）
    paths = {
        '/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'
        '/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'
        '/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'
    };
    for i = 1:length(paths)
        addpath(genpath(paths{i}));
    end

    %% 角度列表
    angle_matrix = -angle_deg:angle_deg;

    %% 创建 mask
    mask_filename = ['origin_' flag '.mat'];
    mask_path = fullfile(FolderPath, subfolders{3}, mask_filename);
    if exist(mask_path, 'file')
        load(mask_path);  % usnr_mask
    else
        usnr_mask = zeros(95,112,117);
        switch flag
            case {'front', 'X'}, usnr_mask(:,56,:) = 1;
            case {'left',  'Z'}, usnr_mask(48,:,:) = 1;
            case {'rot',   'Y'}, usnr_mask(:,:,87) = 1;
            otherwise, error('未知的 flag 值，请使用 "front"、"left"、"rot" 或 "X"、"Y"、"Z"');
        end
        save(mask_path, 'usnr_mask');
    end

    %% 主循环
    for angle = angle_matrix
        subFolder = sprintf('%s_%d', flag, angle);
        subFolderPath = fullfile(FolderPath, subFolder);

        if ~exist(subFolderPath, 'dir')
            fprintf('文件夹不存在，跳过: %s\n', subFolderPath); 
            continue;
        end

        % BASIS_B1m 和 BASIS_LOSS_N 的新路径（mag 文件夹中）
        magFolder = fullfile(FolderPath, subFolder, 'mag');

        try
            load(fullfile(magFolder, 'BASIS_B1m.mat'));
            load(fullfile(magFolder, 'BASIS_LOSS_N.mat'));
        catch ME
            fprintf('加载 BASIS_B1m/BASIS_LOSS_N 时出错: %s\n', ME.message); 
            continue;
        end

        % 加载 object_def.mat
        try
            load(fullfile(subFolderPath, 'object_def.mat'));
        catch ME
            fprintf('加载 object_def.mat 出错: %s\n', ME.message); 
            continue;
        end

        % 查找并加载 material 文件
        matFiles = dir(fullfile(subFolderPath, '*material*.mat'));
        if isempty(matFiles)
            fprintf('未找到 material 文件。\n'); 
            continue;
        end

        try
            load(fullfile(subFolderPath, matFiles(1).name));  % 加载 pan_mask
        catch ME
            fprintf('加载文件 %s 出错: %s\n', matFiles(1).name, ME.message); 
            continue;
        end

        

        usnr_mask_pan = usnr_mask .* pan_mask;
        LOSS = LOSS(1:size(B1m,2), 1:size(B1m,2));

        % 传入计算函数，并传递 abs_rf 和 abs_q 参数
        [usnr_map, b1_map, loss_map, rf_rate] = compute_sos_snr_general(B1m, LOSS, usnr_mask_pan, idxS, abs_rf, abs_q);

        fprintf('当前系统时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

        
        % 保存部分省略，保持原样 ...
        save_items = {
            'usnr_map',         subfolders{1}, 'USNR',     usnr_map;
            'loss_map',         subfolders{2}, 'loss',     loss_map;
            'usnr_mask_pan',    subfolders{3}, 'mask',     usnr_mask_pan;
            'b1_map',           subfolders{4}, 'B1',       b1_map;
            'rf_rate',          subfolders{5}, 'rfrate',   rf_rate;
        };

        for k = 1:size(save_items, 1)
            varname = save_items{k,1};
            subfolder = save_items{k,2};
            prefix = save_items{k,3};
            value = save_items{k,4}; %#ok<NASGU>

            outfile = fullfile(FolderPath, subdir, subfolder, sprintf('sos_%s_%s_%d.mat', prefix, flag, angle));
            save(outfile, varname);
        end

        correction_items = {
            'usnr_map', subfolders{1}, 'USNR', shift_nonzero_region(usnr_map, flag, -angle);
            'loss_map', subfolders{2}, 'loss', shift_nonzero_region(loss_map, flag, -angle);
            'b1_map',   subfolders{4}, 'B1',   shift_nonzero_region(b1_map, flag, -angle);
        };

        for k = 1:size(correction_items, 1)
            varname = correction_items{k,1};
            subfolder = correction_items{k,2};
            prefix = correction_items{k,3};
            value = correction_items{k,4};

            outfile = fullfile(FolderPath, subdir, subfolder, sprintf('sos_%s_%s_%d_correction.mat', prefix, flag, angle));
            save(outfile, 'value');
        end

        clear b1_map loss_map B1m usnr_map usnr_mask_pan
        cd('..');
    end
end

function tf = is_writable(path)
    tf = false;
    testfile = tempname(path);
    fid = fopen(testfile, 'w');
    if fid ~= -1
        tf = true;
        fclose(fid); delete(testfile);
    end
end
