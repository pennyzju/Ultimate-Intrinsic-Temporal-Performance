% 启动交互式切片查看器
V = alpha;
% 假设 V 是你的三维矩阵
if any(isnan(V(:))) || any(isinf(V(:)))
    fprintf('检测到矩阵中含有 NaN 或 Inf，正在修复...\n');
    
    % 将 NaN 替换为 0 (或者替换为矩阵中的最小值)
    V(isnan(V)) = 0;
    
    % 将 Inf (无穷大) 替换为矩阵中的最大有限值
    % 或者直接替换为 0，视具体物理意义而定
    maxVal = max(V(~isinf(V))); 
    if isempty(maxVal), maxVal = 0; end % 防止全是 Inf 的极端情况
    V(isinf(V)) = maxVal;
end
h = orthosliceViewer(V);

% 提示：
fprintf('操作说明：\n');
fprintf('1. 在窗口中点击或拖动十字准星，找到你感兴趣的点。\n');
fprintf('2. 选好后，在命令行输入：pos = h.CrosshairLocation\n');
fprintf('3. pos 的结果即为 [Column, Row, Plane]，对应矩阵的 [Y, X, Z]\n');