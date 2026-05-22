function Basis_Pan_script(tesla, out_dir, flag, pan_step)
    % Basis_Pan_script: 计算给定方向和步长下的ubasis结果
    % 输入参数：
    %   out_dir  - 输出目录，如 '/data2/jiaxinli/pan/20250401_7T_UISNR_output/'
    %   flag     - 方向标志 'X', 'Y', 或 'Z'
    %   pan_step - 扫描角度步长，例如 2
    
        addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
    
        clearvars -except out_dir flag pan_step
        close all;
    
        public_out_dir = fullfile(out_dir, 'public'); % 公开结果目录
        if ~exist(public_out_dir, 'dir') 
            mkdir(public_out_dir); 
        end
    
        filename1 = fullfile(public_out_dir, 'ubasis_options.mat');
        if ~exist(filename1, 'file')
            default_config(public_out_dir,tesla);
        end
    
        filename2 = fullfile(public_out_dir, 'default_material_maps.mat');
        if ~exist(filename2, 'file')
            default_material_maps(public_out_dir);
        end
        load(filename2); % 加载 epsilon_r, sigma_e, rho 等变量
    
        % generate dipoles
        filename3 = fullfile(public_out_dir, 'dist_dipolesObject.mat');
        mask = double( epsilon_r > 1 & sigma_e > 0 & rho > 0 );
        if ~exist(filename3, 'file')
            ubasis_comp_dist_dipoleObject(mask, r, public_out_dir);
        else
            load(filename3); % 加载 mask_dip
            if any(size(mask_dip) ~= size(mask))
                ubasis_comp_dist_dipoleObject(mask, r, public_out_dir);
            end
        end
    
        pan_matrix = -pan_step:1:-1;
    
        for i = 1:length(pan_matrix)
            folderName = sprintf('%s_%d', flag, pan_matrix(i));
            folderPath = fullfile(out_dir, folderName);
    
            if ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end
    
            fileName = sprintf('%s_%d_material_maps.mat', flag, pan_matrix(i));
            pan_dis = pan_matrix(i);
    
            fprintf('PREPARE ...\n');
            pan_maxit_prepare(public_out_dir, folderPath, fileName, flag, pan_dis);
            fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
    
            fprintf('COMPUTE INCIDENT BASIS ...\n');
            ubasis_comp_inc_basis_pan(public_out_dir, folderPath, fileName);
            fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
    
            fprintf('SOLVE ...\n');
            ubasis_solve_pan(public_out_dir, folderPath, fileName);
            fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
    
            fprintf('COMPUTE SCATTERED FIELDS ...\n');
            ubasis_comp_fields_pan(public_out_dir, folderPath, fileName);
            fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
        end
    
        fprintf('%s basis finish ...\n', flag);
    end
    