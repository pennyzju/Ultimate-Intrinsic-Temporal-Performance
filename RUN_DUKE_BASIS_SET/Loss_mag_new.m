function Loss_mag_new(out_dir, flag, angle_matrix)
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
    fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
    close all;

    % COMPUTATION OPTIONS
    public_out_dir = fullfile(out_dir, 'public'); % 公开结果目录
    if ~exist(public_out_dir, 'dir')
        mkdir(public_out_dir);
    end

    %angle_deg = 5;
    %angle_matrix = -angle_deg : 1 : angle_deg;
    %angle_matrix = 3;
    subfolder = 'mag';

    for angle = angle_matrix
        % 在程序开头加入
        poolobj = gcp('nocreate');  % 获取当前并行池对象（如果存在）
        if ~isempty(poolobj)
            delete(poolobj);  % 关闭现有并行池
        end
        folderName = sprintf('%s_%d', flag, angle);
        fprintf('%s is computing ...\n', folderName);

        folderPath = fullfile(out_dir, folderName);
        if ~exist(folderPath, 'dir')
            fprintf('跳过不存在的子文件夹: %s\n', folderPath);
            continue;
        end

        fileName = sprintf('%s_%d_material_maps.mat', flag, angle);

        fprintf('COMPUTE SCATTERED FIELDS ...\n');
        ubasis_comp_fields_new_mag(public_out_dir, folderPath, fileName , subfolder);

        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
        fprintf('-----------------------------\n');
    end

    fprintf('%s basis finish ...\n', flag);
end
