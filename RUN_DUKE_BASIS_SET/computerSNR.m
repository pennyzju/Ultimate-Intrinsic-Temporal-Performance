
filepath = "/data/jiaxinli/projects/coil/results/20240601";

filename = 'object_def.mat';
% 使用 fullfile 构建完整路径并加载文件
fullpath = fullfile(filepath, filename);
load(fullpath);

% 或者，使用字符串连接构建完整路径并加载文件
% fullpath = [filepath '/' filename];
% load(fullpath);

% 或者，使用 sprintf 构建完整路径并加载文件
% fullpath = sprintf('%s/%s', filepath, filename);
% load(fullpath);



compute_usnr_unacc_fast_v3(63,80,85,500,filepath,1e-3,mask,"LOSS_N_500")
