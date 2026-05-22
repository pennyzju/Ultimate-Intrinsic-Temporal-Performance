function snr_process_opt_sta_M(B0, subdir,flag, angle_matrix)
    %B0 = 3, subdir = '20250730staopt',flag = 'front',angle_matrix = [-5:1:5];
% 处理全角度的 USNR/B1/loss 映射，并保存结果
% 参数：
%   flag: 'X', 'Y', 'Z' 或自定义如 'front'/'rot'
%   angle_deg: 最大旋转角度（整数）
%   FolderPath: 工作文件夹路径

%% 添加路径(注意放在最前方)
paths = {
    '/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'
    '/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'
    '/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'
};
for i = 1:length(paths)
    addpath(genpath(paths{i}));
end

%% 构建子文件夹
FolderPath = B0_folderpath(B0);
subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map','rf_map'};
for i = 1:length(subfolders)
    folder_name = fullfile(FolderPath,subdir, subfolders{i});
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
end

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
%% 加载 rf 和 mask 索引
R = load(fullfile(FolderPath,subdir, subfolders{5}, sprintf('rf_%s_0.mat', flag)));
rf_all = R.rf_map;
%mask_rf = load(fullfile(FolderPath, subdir, subfolders{3}, sprintf('mask_%s_0.mat', flag)));
%subFolder = sprintf('%s_0', flag);
mask_path = fullfile(FolderPath, subdir, subfolders{3}, sprintf('rf_mask_cal_%s.mat',flag));
tmp = load(mask_path, 'rf_mask_cal');
rf0mask = tmp.rf_mask_cal.*usnr_mask;
idxS0 =find(rf0mask>0);
s = size(idxS0);
fprintf('idxS0 的大小为 %d 行 × %d 列，总元素个数: %d\n', s(1), s(2), numel(idxS0));


    %% 主计算循环
    for i = 1:length(angle_matrix)
        angle = angle_matrix(i);
        subFolder = sprintf('%s_%d', flag, angle);
        subFolderPath = fullfile(FolderPath, subFolder);
        try
            B = load(fullfile(subFolderPath, 'BASIS_B1m.mat'));
            L = load(fullfile(subFolderPath, 'BASIS_LOSS_N.mat'));
            O = load(fullfile(subFolderPath, 'object_def.mat'));
            M = load(fullfile(subFolderPath, sprintf('%s_%d_material_maps', flag, angle)));
            B1m = B.B1m;
            LOSS = L.LOSS(1:size(B1m,2), 1:size(B1m,2));
            idxS = O.idxS;
            mask = M.rotated_mask;
        
        catch ME
            logmsg(['基础文件加载失败: ' ME.message]);
            continue;
        end

        usnr_mask_rot = usnr_mask .* mask;
    %     fprintf('B1m size: %s\n', mat2str(size(B1m)));
    % fprintf('rf_all size: %s\n', mat2str(size(rf_all)));
    % fprintf('idxS size: %s\n', mat2str(size(idxS)));
    % fprintf('idxS0 size: %s\n', mat2str(size(idxS0)));
    % fprintf('usnr_mask_rot size: %s\n', mat2str(size(usnr_mask_rot)));
    % fprintf('LOSS size: %s\n', mat2str(size(LOSS)));
        LOSS = LOSS(2:2:size(LOSS,1),2:2:size(LOSS,1));
        B1m = B1m(:,2:2:size(B1m,2));

        [b1_map, loss_map, usnr_map] = compute_usnr_rf_sta(B1m, rf_all, idxS, idxS0 , usnr_mask_rot, LOSS);

        % 保存原始 map
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('sta_USNR_%s_%d.mat', flag, angle)), 'usnr_map');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('sta_loss_%s_%d.mat', flag, angle)), 'loss_map');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('sta_B1_%s_%d.mat', flag, angle)), 'b1_map');
        

        % 保存角度校正 map
        [p1, p2] = define_rotation_axis(flag);
        usnr_map_correction= rotate_matrix(usnr_map, p1, p2, -angle);
        loss_map_correction= rotate_matrix(loss_map, p1, p2, -angle);
        b1_map_correction= rotate_matrix(b1_map, p1, p2, -angle);
        save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('sta_USNR_%s_%d_correction.mat', flag, angle)), 'usnr_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('sta_loss_%s_%d_correction.mat', flag, angle)), 'loss_map_correction');
        save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('sta_b1_%s_%d_correction.mat', flag, angle)), 'b1_map_correction');

        logmsg(['完成角度 ' num2str(angle)]);
        %clear b1_map loss_map B1m usnr_map usnr_mask_rot

    end
end
%% 工具函数：日志打印
function logmsg(msg)
    t = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s] %s\n', t, msg);
end
