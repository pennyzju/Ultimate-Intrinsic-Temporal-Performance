
outputfolder = 'H:\UISNR\20240801_UISNR_output\front_1';
data1 = 'BASIS_B1m.mat';
data2 = 'object_def.mat';

load(fullfile(outputfolder, data1));
load(fullfile(outputfolder, data2));


dim = [95, 112, 117];
M = zeros(dim);
M(a) = B1m(:,1);
visualize_3D_slices_correctly(abs(M));

%% % 假设已有复数磁场矩阵 M
V = abs(M);  % 计算磁场强度
V(V < 1e-20) = NaN;
visualize_3D_slices_correctly(V);
% 构建网格
[x, y, z] = meshgrid(1:size(V,2), 1:size(V,1), 1:size(V,3));

% 绘图开始
figure; hold on;

% ========== 1. 外部壳体：isosurface ==========
outer_thresh = 5e-20;
p = patch(isosurface(x, y, z, V, outer_thresh));
isonormals(x, y, z, V, p);
set(p, 'FaceColor', [0.2 0.2 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % 深蓝色半透明

% ========== 2. 内部剖面：slice ==========
xslice = round(size(V,2)/2);
yslice = round(size(V,1)/2);
zslice = round(size(V,3)/2);

hs = slice(x, y, z, V, xslice, yslice, zslice);
set(hs, 'EdgeColor', 'none', 'FaceAlpha', 1);  % slice 不透明，显示清晰
shading interp;

% ========== 3. 美化 ==========
axis tight; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3); camlight headlight; lighting gouraud;
colormap(jet);  % 或 colormap(turbo)
colorbar;
title('|B1| Distribution with Slice View');
% 隐藏坐标轴
axis off;
% ========== 4. 可选导出 ==========
% saveas(gcf, 'magnetic_field_visualization.png');



