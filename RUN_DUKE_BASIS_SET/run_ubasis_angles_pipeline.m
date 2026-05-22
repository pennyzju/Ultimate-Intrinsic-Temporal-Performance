%function run_ubasis_angles_pipeline(B0_str, flag, angle_matrix)
    %RUN_UBASIS_ANGLES Run ubasis pipeline for given field strength, rotation flag and angles
    %
    % INPUTS:
    %   B0_str      - field strength string, e.g., '3T', '1p5T', '7T', '10p5T', '14T'
    %   flag        - 'left', 'rot', or 'front'
    %   angle_matrix- array of rotation angles in degrees, e.g., -2:1:-1
    B0_str = '3T'; 
    flag = 'left';
    angle_matrix = 0;
    %% 0. Setup paths
    clearvars -except B0_str flag angle_matrix;
    close all;
    
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
    
    %% 1. Map B0_str -> numeric value
    B0_map = containers.Map({'1p5T','3T','7T','10p5T','14T'}, {1.5,3,7,10.5,14});
    
    if ~isKey(B0_map, B0_str)
        error('Unsupported B0_str: %s', B0_str);
    end
    B0 = B0_map(B0_str);
    
    %% 2. Define output directories
    out_dir = B0_folderpath_new(B0); % 用户自定义函数，需要返回路径
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    
    public_out_dir = fullfile(out_dir, 'public');
    if ~exist(public_out_dir, 'dir')
        mkdir(public_out_dir);
    end
    
    %% 3. Default options
    opt_file = fullfile(public_out_dir, 'ubasis_options.mat');
    if ~exist(opt_file, 'file')
        default_config(public_out_dir, B0);
    end
    
    %% 4. Default material maps
    mat_file = fullfile(public_out_dir, 'default_material_maps.mat');
    if ~exist(mat_file, 'file')
        default_material_maps(public_out_dir, B0_str);
    end
    load(mat_file);
    
    %% 5. Generate dipoles
    dip_file = fullfile(public_out_dir, 'dist_dipolesObject.mat');
    mask = double(epsilon_r>1 & sigma_e>0 & rho>0);
    
    if ~exist(dip_file,'file')
        ubasis_comp_dist_dipoleObject(mask, r, public_out_dir);
    else
        load(dip_file);
        if ~isequal(size(mask_dip), size(mask))
            ubasis_comp_dist_dipoleObject(mask, r, public_out_dir);
        end
    end
    
   %% 6. Loop over angles
for i = 1:length(angle_matrix)
    angle = angle_matrix(i);

    folderName = sprintf('%s_%d', flag, angle);
    folderPath = fullfile(out_dir, folderName);

    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    fileName = sprintf('%s_%d_material_maps.mat', flag, angle);

    fprintf('\n===== Computing %s =====\n', folderName);

    % 记录循环开始时间
    tStart = tic;

    fprintf('PREPARE ...\n');
    rotate_maxit_prepare(public_out_dir, folderPath, fileName, flag, angle);
    fprintf('Current time: %s\n', datestr(now));

    fprintf('COMPUTE INCIDENT BASIS ...\n');
    ubasis_comp_inc_basis_new(public_out_dir, folderPath, fileName);
    fprintf('Current time: %s\n', datestr(now));

    fprintf('SOLVE ...\n');
    ubasis_solve_new(public_out_dir, folderPath, fileName);
    fprintf('Current time: %s\n', datestr(now));

    fprintf('COMPUTE SCATTERED FIELDS ...\n');
    ubasis_comp_fields_new(public_out_dir, folderPath, fileName);
    fprintf('Current time: %s\n', datestr(now));

    % 计算循环耗时（小时）
    elapsedTime_h = toc(tStart) / 3600;
    fprintf('=== Loop for angle %d finished. Elapsed time: %.3f hours ===\n', angle, elapsedTime_h);

end

    
    fprintf('\n[%s] basis computation finished for angles %s\n', flag, mat2str(angle_matrix));
    
    %end
    
    