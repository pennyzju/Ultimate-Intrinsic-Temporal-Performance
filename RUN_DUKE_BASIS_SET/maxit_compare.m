% 假设有两个矩阵 A 和 B
A = sigma_e;
B = rho;

% 获取两个矩阵的非零元素的索引
[idxA_row, idxA_col] = find(A);
[idxB_row, idxB_col] = find(B);

% 将索引对转换为字符串数组，以便于比较
idxA = strcat(string(idxA_row), ',', string(idxA_col));
idxB = strcat(string(idxB_row), ',', string(idxB_col));

% 检查索引是否完全相同
if isequal(sort(idxA), sort(idxB))
    disp('矩阵 A 和 B 的非零数据下标完全一样。');
else
    disp('矩阵 A 和 B 的非零数据下标不一样。');
end
