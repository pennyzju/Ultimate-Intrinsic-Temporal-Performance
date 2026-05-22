function basis_script(B0,flag,angle_matrix)
    %flag = 'front'; % 'left' or 'rot' or 'front';
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
    %clear;
    close all;

    %  COMPUTATION OPTIONS
    out_dir = B0_folderpath(B0);
    public_out_dir = fullfile(out_dir, 'public'); % 公开结果目录
    if ~exist(public_out_dir, 'dir') 
        mkdir(public_out_dir); 
    end

    filename1 = sprintf('%s/ubasis_options.mat',public_out_dir);
    if ~exist(filename1, 'file')
        % create default configuration options
        default_config(public_out_dir,B0);
    end

    filename2 = sprintf('%s/default_material_maps.mat',public_out_dir);
    if ~exist(filename2, 'file')
        % create default object/material maps
        default_material_maps(public_out_dir);
    end
    load(filename2);

    % % generate dipoles
    filename3 = sprintf('%s/dist_dipolesObject.mat',public_out_dir);
    mask = double( epsilon_r>1 & sigma_e>0 & rho>0 );
    if ~exist(filename3,'file')
    ubasis_comp_dist_dipoleObject(mask,r,public_out_dir);
    else
        load(filename3);
        if size(mask_dip) ~= size(mask)
            ubasis_comp_dist_dipoleObject(mask,r,public_out_dir);
        end
    end

    % 循环生成文件名，并写入数据
    for i = 1:length(angle_matrix)
        % 在程序开头加入
        poolobj = gcp('nocreate');  % 获取当前并行池对象（如果存在）
        if ~isempty(poolobj)
            delete(poolobj);  % 关闭现有并行池
        end
        % 构建文件夹名称
        folderName = sprintf('%s_%d', flag, angle_matrix(i));
        fprintf('正在计算：%s\n', folderName);  

        % 构建完整的文件夹路径
        folderPath = fullfile(out_dir, folderName);

        % 检查文件夹是否存在
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end


        % 生成文件名
        fileName = sprintf('%s_%d_material_maps.mat', flag, angle_matrix(i));
        angle = angle_matrix(i); % 旋转角度(最大可以45度)
        % PREPARE
        fprintf('PREPARE ...\n');
        rotate_maxit_prepare(public_out_dir,folderPath,fileName, flag, angle);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

        %compute incident fields
        fprintf('COMPUTE INCIDENT BASIS ...\n');
        ubasis_comp_inc_basis_new(public_out_dir,folderPath,fileName);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

        % SOLVE FOR INDUCED CURRENTS (J-VIE)
        fprintf('SOLVE ...\n');
        ubasis_solve_new(public_out_dir,folderPath,fileName);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

        % COMPUTE SCATTERED FIELDS
        fprintf('COMPUTE SCATTERED FIELDS ...\n');
        ubasis_comp_fields_new(public_out_dir,folderPath,fileName);
        fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));

    end
    fprintf('%s basis finish ...\n', flag);

end









