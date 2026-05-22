clear
% 获取当前工作目录
%FolderPath  = pwd;
folder_name = './loss';  % 使用 '../' 表示与当前文件夹同级的目录,'./'表示当前文件夹
if ~exist(folder_name, 'dir')
    mkdir(folder_name);  % 如果文件夹不存在，则创建它
end
addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'))


% 获取文件夹下所有子文件夹的信息
subFolders = dir(pwd);

% 遍历所有子文件夹
for i = 1:length(subFolders)
    % 检查是否为文件夹并且名称中包含下划线
    if subFolders(i).isdir && contains(subFolders(i).name, '_')
        % 获取符合条件的子文件夹路径
        subFolderPath = fullfile(pwd, subFolders(i).name);
        fprintf('进入子文件夹: %s\n', subFolderPath);
        
        % 进入子文件夹
        cd(subFolderPath);
        
        % 加载指定的 .mat 文件
        try
            load('BASIS_B1m');
            load('object_def');
            load('BASIS_LOSS_N');
            fprintf('成功加载基础文件。\n');
        catch ME
            fprintf('加载基础文件时出错: %s\n', ME.message);
        end
        
        % 加载文件名中包含 'material' 的 .mat 文件
        matFiles = dir('*material*.mat');  % 查找包含 'material' 的文件
        try
            load(matFiles(j).name);
            fprintf('成功加载文件: %s\n', matFiles(j).name);
        catch ME
            fprintf('加载文件 %s 时出错: %s\n', matFiles(j).name, ME.message);
            
        end
        
%         LOSS=LOSS(1:size(B1m,2),1:size(B1m,2));
% 
%         [usnr_map,b1_map,loss_map,rf_rate]=compute_usnr_unacc_fast_v4(B1m,LOSS,rotated_mask,idxS);
%         fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
%         % 获取当前工作目录
        currentFolderPath = pwd;

        % 使用 fileparts 提取当前文件夹名称
        [~, currentFolderName, ~] = fileparts(currentFolderPath);


        Xout_file =  sprintf('../SNR/USNR_%s_full.mat', currentFolderName);
        save(Xout_file, 'usnr_map', '-v7.3');
        Xout_file2 = ['./loss/','loss_',num2str(seq(i)),'_full.mat'];
        save(Xout_file2, 'loss_map', '-v7.3');
        % 返回到原始目录
        cd(folderPath);
    end
end

