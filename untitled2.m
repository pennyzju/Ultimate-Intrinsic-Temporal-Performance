% 假设 mask 和 mask_dip2 尺寸一致
kk1 = smooth3(mask);         % mask 平滑
kk2 = smooth3(mask_dip2);    % mask_dip2 平滑

figure;

% 绘制 mask（红色，50%透明）
p1 = patch(isosurface(kk1, 0.1));
set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.7);

hold on;

% 绘制 mask_dip2（灰色，50%透明）
%p2 = patch(isosurface(kk2, 0.1));
%set(p2, 'FaceColor', [0.94,0.94,0.94], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% 视觉效果
daspect([1 1 1]);  % 坐标比例一致
%view(3);            % 3D 视图
%view([1,1,1]);
view([110  18]);
%view([90    0]);%正视图
%view([90    90]);%俯视图
%view([180    0]);%侧视图
camlight;           % 加光源
lighting gouraud;   % 光照效果
axis off;           % 不显示坐标轴
% 保存为 PNG
filename = 'head_surface.png';
%saveas(gcf, filename);   % 方法1：saveas
% 或者使用 exportgraphics（更高质量）
exportgraphics(gcf, filename, 'Resolution',300);