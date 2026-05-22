% 定义源文件完整路径
sourceFile = 'F:\20250201_1_5T_UISNR_output\20250730staopt\SNR_plain\sta_USNR_front_0.mat';

% 定义目标文件完整路径（包含新文件名）
targetFile = 'F:\20250201_1_5T_UISNR_output\20251022seq\SNR_plain\USNR_front_0_2500.mat';

% 检查源文件是否存在
if exist(sourceFile, 'file')
    % 执行复制并重命名
    [status, msg] = copyfile(sourceFile, targetFile);
    
    if status
        fprintf('成功：文件已复制并重命名为 USNR_left_0_2500.mat\n');
    else
        fprintf('失败：错误信息 - %s\n', msg);
    end
else
    fprintf('错误：源文件不存在，请检查路径是否正确。\n');
end

%% 
% 定义源文件完整路径
sourceFile = 'F:\20250201_1_5T_UISNR_output\20250730staopt\B1_map\sta_B1_left_0.mat';

% 定义目标文件完整路径（包含新文件名）
targetFile = 'F:\20250201_1_5T_UISNR_output\20251022seq\B1_map\B1_left_0_2500.mat';

% 检查源文件是否存在
if exist(sourceFile, 'file')
    % 执行复制并重命名
    [status, msg] = copyfile(sourceFile, targetFile);
    
    if status
        fprintf('成功：文件已复制并重命名为 USNR_left_0_2500.mat\n');
    else
        fprintf('失败：错误信息 - %s\n', msg);
    end
else
    fprintf('错误：源文件不存在，请检查路径是否正确。\n');
end
%% 
% 定义源文件完整路径
sourceFile = 'H:\UISNR\20240801_UISNR_output\20250730staopt\SNR_plain\sta_USNR_front_0.mat';

% 定义目标文件完整路径（包含新文件名）
targetFile = 'F:\20250201_1_5T_UISNR_output\20251022seq\SNR_plain\USNR_front_0_2500.mat';

% 检查源文件是否存在
if exist(sourceFile, 'file')
    % 执行复制并重命名
    [status, msg] = copyfile(sourceFile, targetFile);
    
    if status
        fprintf('成功：文件已复制并重命名为 USNR_left_0_2500.mat\n');
    else
        fprintf('失败：错误信息 - %s\n', msg);
    end
else
    fprintf('错误：源文件不存在，请检查路径是否正确。\n');
end

