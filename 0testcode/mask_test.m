clear

folder_name = './mask_test';  % 使用 '../' 表示与当前文件夹同级的目录,'./'表示当前文件夹
if ~exist(folder_name, 'dir')
    mkdir(folder_name);  % 如果文件夹不存在，则创建它
end
addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'))
addpath(genpath('H:\UISNR\20240801_UISNR_output'));
flag = 'front'; % 'left' or 'rot' or 'front';
angle_deg =5;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
angle_matrix = -angle_deg:1:angle_deg;
% 删除角度矩阵中的0
angle_matrix(angle_matrix == 0) = [];  % 使用逻辑索引去除0

for i = 1:length(angle_matrix)
    % 构建文件夹名称
    subFolder = sprintf('%s_%d', flag, angle_matrix(i));
    subFolderPath = fullfile(pwd, subFolder);
    fprintf('进入子文件夹: %s\n', subFolderPath);
     % 进入子文件夹
     cd(subFolderPath);
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
            end
        end
        % 创建一个与 A 同尺寸的零矩阵mask
        mask_zeroed = zeros(m,n,p);
        % 只保留中间层，将其他层数据置零
        mask_zeroed(ceil(m / 2),:,:) = 1;
        mask_zeroed(:,ceil(n / 2),:) = 1;
        mask_zeroed(:,:,ceil(p / 2)) = 1;
        % 旋转 mask_zeroed
        mask_zeroed = rotate_matrix(mask_zeroed, point1, point2, angle_matrix(i));

        % 只保留中间层，将其他层数据置零
        mask_zeroed = mask_zeroed.*rotated_mask;
        
        save(sprintf('../mask_test/%s_%d.mat', flag, angle_matrix(i)), 'mask_zeroed');     
        % 返回到原始目录
        cd('..');
end




