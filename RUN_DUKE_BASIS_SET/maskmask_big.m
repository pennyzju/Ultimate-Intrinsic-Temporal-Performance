% 设定文件目录
public_out_dir = 'F:\20240728test\output';

% 获取目录中所有包含 material_maps 的 .mat 文件
files = dir(fullfile(public_out_dir, '*material_maps.mat'));

% 初始化累加的 mask 矩阵
mask_sum = [];

% 循环遍历每个文件
for i = 1:length(files)
    % 获取文件的完整路径
    filePath = fullfile(files(i).folder, files(i).name);
    
    % 加载文件
    data = load(filePath);
    
    % 检查文件中是否包含 mask 矩阵
    if isfield(data, 'rotated_mask')
        % 如果 mask_sum 为空，则初始化为相同大小的零矩阵
        if isempty(mask_sum)
            mask_sum = zeros(size(data.mask));
        end
        
        % 累加 mask 矩阵
        mask_sum = mask_sum + data.rotated_mask;
    else
        fprintf('文件 %s 中不包含 mask 矩阵。\n', files(i).name);
    end
end

disp('所有 mask 矩阵累加完毕。');

