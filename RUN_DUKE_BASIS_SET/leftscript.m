addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));

clear;
close all;

%  COMPUTATION OPTIONS
public_out_dir = '/data/jiaxinli/projects/20240801_UISNR_output/public';%结果输出目录

filename1 = sprintf('%s/ubasis_options.mat',public_out_dir);
if ~exist(fullfile(public_out_dir, 'ubasis_options.mat'), 'file')
% create default configuration options
    default_config(public_out_dir);
end

filename2 = sprintf('%s/default_material_maps.mat',public_out_dir);
% create default object/material maps
if ~exist(fullfile(public_out_dir, 'default_material_maps.mat'), 'file')
    default_material_maps(public_out_dir);
end
load(filename2);

out_dir = '/data/jiaxinli/projects/20240801_UISNR_output/';%结果输出目录
% 检查目录是否存在，如果不存在则创建它
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% % generate dipoles
filename3 = sprintf('%s/dist_dipolesObject.mat',public_out_dir);
mask = double( epsilon_r>1 & sigma_e>0 & rho>0 );
if ~exist(filename3,'file')
   ubasis_comp_dist_dipoleObject(mask,r,public_out_dir);
else
    load(filename3);
    if size(mask_dip) ~= size(mask)
        ubasis_comp_dist_dipoleObject(mask,r,public_out_dir);
    end
end

flag = 'left'; % 'left' or 'rot' or 'front';
angle_deg =5 ;
% 生成从 -angle 到 angle，间隔为 1 度的角度矩阵
angle_matrix = 1:1:angle_deg;


%{
 filename4 = sprintf('%s/BASIS_Einc.mat',public_out_dir);
filename5 = sprintf('%s/BASIS_Hinc.mat',public_out_dir);
% 检查两个文件是否都存在
if exist(filename4, 'file') == 2 && exist(filename5, 'file') == 2
    fprintf('INCIDENT BASIS EXIST \n');
else
    % COMPUTE INCIDENT BASIS
    fprintf('COMPUTE INCIDENT BASIS ...\n');
    ubasis_comp_inc_basis(public_out_dir);
end 
%}



% 循环生成文件名，并写入数据
for i = 1:length(angle_matrix)
    % 构建文件夹名称
    folderName = sprintf('%s_%d', flag, angle_matrix(i));

    % 构建完整的文件夹路径
    folderPath = fullfile(out_dir, folderName);

    % 检查文件夹是否存在
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    
    % 生成文件名
    fileName = sprintf('%s_%d_material_maps.mat', flag, angle_matrix(i));
    angle = angle_matrix(i); % 旋转角度(最大可以45度)
    % PREPARE
    fprintf('PREPARE ...\n');
    rotate_maxit_prepare(public_out_dir,folderPath,fileName, flag, angle);
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

    %ubasis_prepare(public_out_dir,folderPath, fileName, flag);
    %compute incident fields
    fprintf('COMPUTE INCIDENT BASIS ...\n');
    ubasis_comp_inc_basis_new(public_out_dir,folderPath,fileName);
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

    % SOLVE FOR INDUCED CURRENTS (J-VIE)
    fprintf('SOLVE ...\n');
    ubasis_solve_new(public_out_dir,folderPath,fileName);
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

    % COMPUTE SCATTERED FIELDS
    fprintf('COMPUTE SCATTERED FIELDS ...\n');
    ubasis_comp_fields_new(public_out_dir,folderPath,fileName);
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

end
fprintf('%s basis finish ...\n', flag);










