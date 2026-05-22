clear
%% 初始化文件夹
%FolderPath = '/data2/jiaxinli/pan/20250401_7T_UISNR_output';
FolderPath = '/data/jiaxinli/projects/20240801_UISNR_output/';
subfolders = {'SNR_plain', 'loss_plain', 'mask_plain'};

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

flag = 'left'; % 'left' or 'rot' or 'front';
angle_deg =5;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
%angle_matrix = 1:1:angle_deg;
angle_matrix = -angle_deg:1:-1;
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

    % 根据 flag 设置对应层面为1
    switch flag
        case 'front'
            usnr_mask(:,56,:) = 1;     % Y = 56
        case 'left'
            usnr_mask(48,:,:) = 1;     % X = 48
        case 'rot'
            usnr_mask(:,:,87) = 1;     % Z = 87
        otherwise
            error('未知的 flag 值，请使用 "front"、"left" 或 "rot"');
    end

    % 保存 mask
    save(mask_path, 'usnr_mask');
    fprintf('已生成并保存 usnr_mask：%s\n', mask_path);
end

%计算USNR

for i = 1:length(angle_matrix)
    % 构建文件夹名称
    subFolder = sprintf('%s_%d', flag, angle_matrix(i));
    subFolderPath = fullfile(FolderPath, subFolder);
    % 判断文件夹是否存在
    if exist(subFolderPath, 'dir')
        cd(subFolderPath);  % 进入文件夹
        fprintf('进入子文件夹: %s\n', subFolderPath);
    else
        fprintf('文件夹不存在，跳过: %s\n', subFolderPath);
        continue;  % 跳过本轮循环
    end
    
    %获取数据维度
    m = 95;
    n = 112;
    p = 117;

    % 根据flag的值选择旋转轴
    if strcmp(flag, 'left')
        point1 = [0, n/2, 30];
        point2 = [m, n/2, 30];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    elseif strcmp(flag, 'front')
        point1 = [m/2, 0, 30];
        point2 = [m/2, n, 30];
    elseif strcmp(flag, 'rot')
        point1 = [m/2, n/2, 0];
        point2 = [m/2, n/2, p];
    else
        error('Invalid flag value. Flag should be ''left'', ''front'', or ''rot''.');
    end
    
        % 加载指定的 .mat 文件
        try
            load('BASIS_B1m');
            load('object_def');
            load('BASIS_LOSS_N');
            fprintf('成功加载基础文件。\n');
        catch ME
            fprintf('加载基础文件时出错: %s\n', ME.message);
            continue;  % 跳过本轮循环
        end
        
        % 加载文件名中包含 'material' 的 .mat 文件
        % 查找当前目录中所有包含 'material' 的 .mat 文件
        matFiles = dir('*material*.mat');

        % 检查 matFiles 是否为空
        if isempty(matFiles)
            disp('没有找到包含 "material" 的 .mat 文件。');
        else
            try
                % 加载第一个找到的文件
                load(matFiles(1).name);
                fprintf('成功加载文件: %s\n', matFiles(1).name);
            catch ME
                fprintf('加载文件 %s 时出错: %s\n', matFiles(1).name, ME.message);
                continue;  % 跳过本轮循环
            end
        end

        % 与原始矩阵求交集
        usnr_mask_rot = usnr_mask.*rotated_mask;       

        LOSS=LOSS(1:size(B1m,2),1:size(B1m,2));

        [usnr_map,b1_map,loss_map]=compute_usnr_unacc_fast_v4(B1m,LOSS,usnr_mask_rot,idxS);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

        % 保存 USNR map
        Xout_file1 = fullfile(FolderPath, subfolders{1}, sprintf('USNR_%s_%d.mat', flag, angle_matrix(i)));
        save(Xout_file1, 'usnr_map');  % 默认保存格式

        % 保存 loss map
        Xout_file2 = fullfile(FolderPath, subfolders{2}, sprintf('loss_%s_%d.mat', flag, angle_matrix(i)));
        save(Xout_file2, 'loss_map');  % 默认保存格式

        % 保存 mask
        Xout_file3 = fullfile(FolderPath, subfolders{3}, sprintf('mask_%s_%d.mat', flag, angle_matrix(i)));
        save(Xout_file3, 'usnr_mask_rot');  % 默认保存格式

        usnr_map_correction= rotate_matrix(usnr_map, point1, point2, -angle_matrix(i));
        Xout_file4 = fullfile(FolderPath, subfolders{1}, sprintf('USNR_%s_%d_correction.mat', flag, angle_matrix(i)));
        save(Xout_file4, 'usnr_map_correction');
        
        loss_map_correction= rotate_matrix(loss_map, point1, point2, -angle_matrix(i));
        Xout_file5 = fullfile(FolderPath, subfolders{2}, sprintf('loss_%s_%d_correction.mat', flag, angle_matrix(i)));
        save(Xout_file4, 'loss_map_correction');
        clear usnr_map b1_map loss_map B1m usnr_map_correction loss_map_correction usnr_mask_rot
        
        % 返回到原始目录
        cd('..');
end




