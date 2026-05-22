clear
%% 初始化文件夹
FolderPath = '/data2/jiaxinli/pan/20250416_1_5T_UISNR_output';
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
angle_deg =4;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
%angle_matrix = 2:1:angle_deg;
angle_matrix = 0;
% 删除角度矩阵中的0
%angle_matrix(angle_matrix == 0) = [];  % 使用逻辑索引去除0


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
        usnr_mask_pan = usnr_mask.*pan_mask;       

        LOSS=LOSS(1:size(B1m,2),1:size(B1m,2));

        [usnr_map, b1_map, loss_map, rf_rate, rf_all]=compute_usnr_unacc_fast_v5(B1m,LOSS,usnr_mask_pan,idxS);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
        %subfolders = {'SNR_plain', 'loss_plain', 'mask_plain','B1_map','rf_map'};
        % 定义保存信息：变量名、子文件夹、文件名前缀、变量本体
        save_items = {
            'usnr_map',         subfolders{1}, 'USNR',             usnr_map;
            'loss_map',         subfolders{2}, 'loss',             loss_map;
            'usnr_mask_pan',    subfolders{3}, 'mask',             usnr_mask_pan;
            'b1_map',           subfolders{4}, 'B1',               b1_map;
            'rf_all',           subfolders{5}, 'rf',               rf_all;
            'rf_rate',          subfolders{5}, 'rfrate',           rf_rate;
        };

        for k = 1:size(save_items, 1)
            varname   = save_items{k, 1};
            subfolder = save_items{k, 2};
            prefix    = save_items{k, 3};
            value     = save_items{k, 4};
            
            outfile = fullfile(FolderPath, subfolder, sprintf('%s_%s_%d.mat', prefix, flag, angle_matrix(i)));
            save(outfile, varname);  % 使用变量名字符串保存
        end

        % 保存 correction 后缀的 map
        correction_items = {
            'usnr_map', subfolders{1}, 'USNR', shift_nonzero_region(usnr_map, flag, -angle_matrix(i));
            'loss_map', subfolders{2}, 'loss', shift_nonzero_region(loss_map, flag, -angle_matrix(i));
            'b1_map',   subfolders{4}, 'b1',   shift_nonzero_region(b1_map, flag, -angle_matrix(i));
        };

        for k = 1:size(correction_items, 1)
            varname   = correction_items{k, 1};
            subfolder = correction_items{k, 2};
            prefix    = correction_items{k, 3};
            value     = correction_items{k, 4};
            
            outfile = fullfile(FolderPath, subfolder, sprintf('%s_%s_%d_correction.mat', prefix, flag, angle_matrix(i)));
            save(outfile, 'value');  % 保存为变量名 value（也可以换名）
        end

        % 可选清除部分变量
        clear b1_map loss_map B1m loss_map_anti usnr_map_anti usnr_mask_pan

        
        % 返回到原始目录
        cd('..');
end





