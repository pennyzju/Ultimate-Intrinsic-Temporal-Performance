addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));

clear;
close all;

%  COMPUTATION OPTIONS
out_dir = '/data2/jiaxinli/pan/20250401_7T_UISNR_output/'; % 结果输出目录
public_out_dir = fullfile(out_dir, 'public'); % 公开结果目录
if ~exist(public_out_dir, 'dir') 
    mkdir(public_out_dir); 
end

filename1 = sprintf('%s/ubasis_options.mat',public_out_dir);
if ~exist(fullfile(public_out_dir, 'ubasis_options.mat'), 'file')
% create default configuration options
    default_config_7(public_out_dir);
end

filename2 = sprintf('%s/default_material_maps.mat',public_out_dir);
% create default object/material maps
if ~exist(fullfile(public_out_dir, 'default_material_maps.mat'), 'file')
    default_material_maps(public_out_dir);
end
load(filename2);

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

flag = 'X'; % 'X' or 'Y' or 'Z';
pan_step =2;
pan_matrix = -pan_step : 1 : pan_step;

% 循环生成文件名，并写入数据
for i = 1:length(pan_matrix)
    % 构建文件夹名称
    folderName = sprintf('%s_%d', flag, pan_matrix(i));

    % 构建完整的文件夹路径
    folderPath = fullfile(out_dir, folderName);

    % 检查文件夹是否存在
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    % 生成文件名
    fileName = sprintf('%s_%d_material_maps.mat', flag, pan_matrix(i));
    pan_dis = pan_matrix(i); % 旋转角度(最大可以45度)

    subfolder = 'mag';
    % COMPUTE SCATTERED FIELDS
    fprintf('COMPUTE SCATTERED FIELDS ...\n');
    ubasis_comp_fields_pan_mag(public_out_dir, folderPath, fileName,subfolder);
    
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));


end
fprintf('%s basis finish ...\n', flag);











