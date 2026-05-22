clear
%% 初始化文件夹
FolderPath = '/data2/jiaxinli/pan/20250416_3T_UISNR_output';
%FolderPath = '/data/jiaxinli/projects/20240801_UISNR_output/';
subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map','rf_map'};

for i = 1:length(subfolders)
    folder_name = fullfile(FolderPath, subfolders{i});
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
end

%% 添加路径
paths = {
    '/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'
    '/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'
    '/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'
};

for i = 1:length(paths)
    addpath(genpath(paths{i}));
end

flag = 'Y'; % 'X' or 'Y' or 'Z';
angle_deg =2;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
%angle_matrix = 2:1:angle_deg;
angle_matrix = -angle_deg:1:angle_deg;
% 删除角度矩阵中的0
angle_matrix(angle_matrix == 0) = [];  % 使用逻辑索引去除0


% 创建 mask 
% 生成保存路径（文件名根据 flag 命名）
mask_filename = ['origin_' flag '.mat'];
mask_path = fullfile(FolderPath, subfolders{3}, mask_filename);

% 如果文件已存在，跳过生成
if exist(mask_path, 'file')
    fprintf('文件已存在，直接加载：%s\n', mask_path);
    load(mask_path);
else
    % 创建空的 mask 矩阵
    usnr_mask = zeros(95,112,117);

    switch flag
        case {'front', 'X'}
            usnr_mask(:,56,:) = 1;     % Y = 56
        case {'left', 'Z'}
            usnr_mask(48,:,:) = 1;     % X = 48
        case {'rot', 'Y'}
            usnr_mask(:,:,87) = 1;     % Z = 87
        otherwise
            error('未知的 flag 值，请使用 "front"、"left"、"rot"或者“X”、“Y”、“Z”');
    end
    

    % 保存 mask
    save(mask_path, 'usnr_mask');
    fprintf('已生成并保存 usnr_mask：%s\n', mask_path);
end

R = load(fullfile(FolderPath,subfolders{5}, sprintf('rf_%s_0.mat', flag)));
rf_all = R.rf_all;
mask_rf = load(fullfile(FolderPath,subfolders{3}, sprintf('mask_%s_0.mat', flag)));
idxS0 = find(mask_rf.usnr_mask_pan);
%% 加载或计算usnr_mask_all
mask_all_path = fullfile(FolderPath, subfolders{3}, 'usnr_mask_all.mat');

if exist(mask_all_path, 'file')
    load(mask_all_path, 'usnr_mask_all');
    logmsg(['已加载交集掩码: ' mask_all_path]);
else
    logmsg('未找到交集掩码，开始计算...');
    usnr_mask_all = usnr_mask;

    for i = 1:length(angle_matrix)
        subFolder = sprintf('%s_%d', flag, angle_matrix(i));
        subFolderPath = fullfile(FolderPath, subFolder);

        if ~isfolder(subFolderPath)
            logmsg(['文件夹不存在，跳过: ' subFolderPath]);
            continue;
        end

        matFiles = dir(fullfile(subFolderPath, '*material*.mat'));
        if isempty(matFiles)
            logmsg(['未找到 material 文件于: ' subFolderPath]);
            continue;
        end

        try
            S = load(fullfile(subFolderPath, matFiles(1).name));
            if isfield(S, 'pan_mask')
                pan_mask = logical(S.pan_mask);
                usnr_mask_pan = usnr_mask.*pan_mask; 
                usnr_mask_all = usnr_mask_all & usnr_mask_pan;
            else
                logmsg(['未找到 pan_mask: ' matFiles(1).name]);
            end
        catch ME
            logmsg(['加载失败: ' ME.message]);
            continue;
        end
    end

    save(mask_all_path, 'usnr_mask_all');
    logmsg(['交集掩码已保存: ' mask_all_path]);
end

%% 主计算循环

for i = 1:length(angle_matrix)
    angle = angle_matrix(i);
    subFolder = sprintf('%s_%d', flag, angle);
    subFolderPath = fullfile(FolderPath, subFolder);
    rfFolderPath = fullfile(FolderPath, 'rf_map');

    try
        B = load(fullfile(subFolderPath, 'BASIS_B1m.mat'));
        L = load(fullfile(subFolderPath, 'BASIS_LOSS_N.mat'));
        O = load(fullfile(subFolderPath, 'object_def.mat'));
        B1m = B.B1m;
        LOSS = L.LOSS(1:size(B1m,2), 1:size(B1m,2));
        idxS = O.idxS;

    catch ME
        logmsg(['基础文件加载失败: ' ME.message]);
        continue;
    end

    gpuInfo = gpuDevice();
    batchSize = min(1024, floor(gpuInfo.AvailableMemory / 8e3));
    [b1_map, loss_map, usnr_map] = compute_usnr_rf_sta(B1m, rf_all, idxS, idxS0 , usnr_mask_all, LOSS);

    % 保存结果
    save_items = {
        'usnr_map',     subfolders{1}, 'USNR';
        'loss_map',     subfolders{2}, 'loss';
        'b1_map',       subfolders{4}, 'B1';
    };

    for k = 1:size(save_items,1)
        v = save_items{k,1}; folder = save_items{k,2}; prefix = save_items{k,3};
        save(fullfile(FolderPath, folder, sprintf('%s_%s_%d.mat', prefix, flag, angle)), v);
    end

    % 保存 correction 后缀
    corr_items = {
        'usnr_map', subfolders{1}, 'USNR';
        'loss_map', subfolders{2}, 'loss';
        'b1_map',   subfolders{4}, 'b1';
    };
    for k = 1:size(corr_items,1)
        vname = corr_items{k,1}; folder = corr_items{k,2}; prefix = corr_items{k,3};
        shifted = shift_nonzero_region(eval(vname), flag, -angle);
        save(fullfile(FolderPath, folder, sprintf('%s_%s_%d_correction.mat', prefix, flag, angle)), 'shifted');
    end

    logmsg(['完成角度 ' num2str(angle)]);
    clear b1_map loss_map usnr_map
end

%% 工具函数：日志打印
function logmsg(msg)
    t = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s] %s\n', t, msg);
end





