% 原始mask大小
mask = zeros(95,112,117);

% 创建与mask同大小的X、Y、Z矩阵
X = zeros(size(mask));
Y = zeros(size(mask));
Z = zeros(size(mask));

% 设置特定层面为1
X(48,:,:) = 1;     % X方向（sagittal面）
Y(:,56,:) = 1;     % Y方向（coronal面）
Z(:,:,87) = 1;     % Z方向（axial面）

X1=epsilon_r.*X;


% 保存为 .mat 文件
save('X.mat', 'X');
save('Y.mat', 'Y');
save('Z.mat', 'Z');

X1=epsilon_r.*X;

pan_epsilon_r = shift_gtone_region(epsilon_r, 'Y', 10);
imagesc(squeeze(pan_epsilon_r(:,:,87)));
%X方向移动，通过Y截面可以观察到；Y方向移动可以通过Z截面观察到，Z方向移动可以通过X截面观察到；
figure;
imagesc(squeeze(epsilon_r(:,:,87)));