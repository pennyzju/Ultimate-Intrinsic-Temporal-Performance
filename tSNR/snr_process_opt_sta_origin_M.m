function snr_process_opt_sta_origin_M(B0, subdir,flag)
% 处理初始位置的 USNR/B1/loss 映射，并保存结果
% 参数：
% flag: 'X', 'Y', 'Z' 或 'front'/'rot'
% angle_deg: 最大旋转角度（整数）
% FolderPath: 工作文件夹路径

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

%% 加载 origin mask
%mask_filename = ['origin_' flag '.mat'];
%subFolder = sprintf('left_0', );
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

rf_mask_cal = usnr_mask .* maskall;

% %% 加载 rf 和 mask 索引
% R = load(fullfile(FolderPath,subdir, subfolders{5}, sprintf('rf_%s_0.mat', flag)));
% rf_all = R.rf_map;
% mask_rf = load(fullfile(FolderPath, subdir, subfolders{3}, sprintf('mask_%s_0.mat', flag)));
% idxS0 = find(mask_rf.usnr_mask_rot);
subFolderPath = fullfile(FolderPath,'left_0','bigmask');
try
    B = load(fullfile(subFolderPath, 'BASIS_B1m.mat'));
    L = load(fullfile(subFolderPath, 'BASIS_LOSS_N.mat'));
    B1m = B.B1m;
    LOSS = L.LOSS(1:size(B1m,2), 1:size(B1m,2));
    
catch ME
    logmsg(['基础文件加载失败: ' ME.message]);
end
LOSS = LOSS(2:2:size(LOSS,1),2:2:size(LOSS,1));
B1m = B1m(:,2:2:size(B1m,2));
[usnr_map, b1_map, loss_map, rf_rate, rf_map] = compute_usnr_unacc_fast_v5(B1m, LOSS, rf_mask_cal, idxS);

% 保存原始 map
save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('sta_USNR_%s_0.mat', flag)), 'usnr_map');
save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('sta_loss_%s_0.mat', flag)), 'loss_map');
save(fullfile(FolderPath,subdir, subfolders{3}, sprintf('rf_mask_cal_%s.mat',flag)), 'rf_mask_cal');
save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('sta_B1_%s_0.mat', flag)), 'b1_map');
save(fullfile(FolderPath,subdir, subfolders{5}, sprintf('rf_%s_0.mat', flag)), 'rf_map');
save(fullfile(FolderPath,subdir, subfolders{5}, sprintf('rfrate_%s_0.mat', flag)), 'rf_rate');
    
% 保存角度校正 map
[p1, p2] = define_rotation_axis(flag);
usnr_map_correction= rotate_matrix(usnr_map, p1, p2, 0);
loss_map_correction= rotate_matrix(loss_map, p1, p2, 0);
b1_map_correction= rotate_matrix(b1_map, p1, p2, 0);
save(fullfile(FolderPath,subdir, subfolders{1}, sprintf('sta_USNR_%s_0_correction.mat', flag)), 'usnr_map_correction');
save(fullfile(FolderPath,subdir, subfolders{2}, sprintf('sta_loss_%s_0_correction.mat', flag)), 'loss_map_correction');
save(fullfile(FolderPath,subdir, subfolders{4}, sprintf('sta_b1_%s_0_correction.mat', flag)), 'b1_map_correction');

logmsg(['完成角度0 ' ]);
clear b1_map loss_map B1m usnr_map usnr_mask_rot

end
%% 工具函数：日志打印
function logmsg(msg)
    t = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s] %s\n', t, msg);
end
