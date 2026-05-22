addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test'));

close all;
out_dir = '/data2/jiaxinli/pan/20250401_7T_UISNR_output/'; % 结果输出目录
public_out_dir = fullfile(out_dir, 'public'); % 公开结果目录


filename2 = sprintf('%s/default_material_maps.mat',public_out_dir);
% create default object/material maps
if ~exist(fullfile(public_out_dir, 'default_material_maps.mat'), 'file')
    default_material_maps(public_out_dir);
end
load(filename2);
checkfolder= fullfile(out_dir, 'B1map');
% 检查文件夹是否存在
if ~exist(checkfolder, 'dir')
    mkdir(checkfolder);
end   

flag = 'X'; % 'X' or 'Y' or 'Z';
pan_step =5;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
%pan_matrix = pan_step:1:-pan_step;
pan_matrix = -pan_step:1:pan_step;


% 循环生成文件名，并写入数据
for i = 1:length(pan_matrix)
    % 构建文件夹名称
    folderName = sprintf('%s_%d', flag, pan_matrix(i));
    fprintf('%s\n', folderName);
    % 构建完整的文件夹路径
    folderPath = fullfile(out_dir, folderName);
    fprintf('%s\n', folderPath);
    % 检查文件夹是否存在
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    load BASIS_B1m;
    load object_def.mat;
    prev_idxS = [];  % 用于保存前一次的 idxS
    % 检查 idxS 是否一致
    if i == 1
        prev_idxS = idxS;
    else
        if ~isequal(idxS, prev_idxS)
            warning('idxS 不一致：第 %d 个文件夹 "%s"', i, folderName);
        end
    end
    dim = [95, 112, 117];
    M = zeros(dim);
    M(idxS) = B1m(:,1);
    visualize_3D_slices_title(abs(M),folderName,checkfolder);

end











