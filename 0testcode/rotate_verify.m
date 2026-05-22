load mask;
flag = 'left'; % 'left' or 'rot' or 'front';
angle_deg =5;

[m, n, p] = size(mask);

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

rotated_mask = rotate_matrix(mask, point1, point2, angle_deg);
mask2 = rotate_matrix(rotated_mask, point1, point2, -angle_deg);

% 计算逐元素差值矩阵
diff_matrix = mask2 - mask;
slice_index = 56;  % 选择切片索引

figure;
% 可视化差值矩阵某一切片
subplot(2,2,1); imagesc(squeeze(mask(ceil(end/2),:,:))); axis image off; title('mask');
subplot(2,2,2); imagesc(squeeze(rotated_mask(ceil(end/2),:,:))); axis image off; title('rotated_mask');
subplot(2,2,3); imagesc(squeeze(mask2(ceil(end/2),:,:))); axis image off; title('mask2');
subplot(2,2,4); imagesc(diff_matrix(:, :, slice_index));  % 显示差值矩阵切片
colorbar;
title('Difference Matrix Slice');

