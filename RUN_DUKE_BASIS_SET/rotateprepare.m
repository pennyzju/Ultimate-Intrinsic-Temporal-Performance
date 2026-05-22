% volumeViewer 3维查看器
% 定义文件的完整路径
addpath(genpath('D:\icloud\iCloudDrive\Project\my-repo'));
file_path = 'D:\output\material_maps';
% 加载 .mat 文件
data = load(file_path);
load('D:\output\dist_dipolesObject.mat')
% 获取 dens_3mm 矩阵
rho = data.rho;
% 获取矩阵的大小
[m, n, p] = size(rho);

% 定义旋转参数
% 前后摆头
% point1 = [0, n/2, 30]; % 旋转轴上的第一个点
% point2 = [m, n/2, 30]; % 旋转轴上的第二个点
% angle_deg =15; % 旋转角度(最大可以15度）
 
% % 左右摆头
% point1 = [m/2, 0, 30]; % 旋转轴上的第一个点
% point2 = [m/2, n, 30]; % 旋转轴上的第二个点
% angle_deg = 12; % 旋转角度(最大可以12度）

%左右旋转
point1 = [m/2, n/2, 0]; % 旋转轴上的第一个点
point2 = [m/2, n/2, p]; % 旋转轴上的第二个点
angle_deg =90; % 旋转角度(最大可以45度）




% 保存旋转后的矩阵到新的 .mat 文件
output_path = 'rotated_DUKE_model_3mm_FLOOD_7T.mat';
save(output_path, 'rotated_rho');

rho1 = rho+2000*mask_dip2;
rotated_rho1 = rotated_rho+2000*mask_dip2;
% 可视化原始矩阵和旋转后的矩阵
figure;

% 可视化原始矩阵中的某一层 (例如中间层)
subplot(2, 2, 1);
imagesc(squeeze(rho1(:, :, round(p / 2))));
title('原始矩阵的中间层');
axis equal;
colorbar;

% 可视化旋转后矩阵中的同一层
subplot(2, 2, 2);
imagesc(squeeze(rotated_rho1(:, :, round(p / 2))));
title('旋转后矩阵的中间层');
axis equal;
colorbar;

% 使用 slice 可视化3D数据
subplot(2, 2, 3);
slice(double(rho1), round(m / 2), round(n / 2), round(p / 2));
title('原始矩阵的3D视图');
shading interp;
colorbar;

subplot(2, 2, 4);
slice(double(rotated_rho1), round(m / 2), round(n / 2), round(p / 2));
title('旋转后矩阵的3D视图');
shading interp;
colorbar;

% 显示图像
drawnow;
