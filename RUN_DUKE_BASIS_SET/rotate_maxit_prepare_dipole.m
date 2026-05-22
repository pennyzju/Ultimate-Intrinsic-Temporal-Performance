function rotate_maxit_prepare_dipole(public_out_dir, out_dir,fileName, flag, angle_deg)
    try
        % 生成文件路径
        filename0 = fullfile(public_out_dir, 'default_material_maps.mat');

        % 检查并创建默认物质映射文件
        if ~exist(filename0, 'file')
            default_material_maps_dipole(public_out_dir);
        end

        % 加载默认物质映射数据
        load(filename0);

        % 获取数据维度
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

        % 旋转矩阵
        rotated_mask_dip2 = rotate_matrix(mask_dip2, point1, point2, angle_deg);
       


        % 将数据写入文件
        save(fullfile(out_dir, fileName), 'dx','rotated_mask_dip2','r');

    catch ME
        % 捕获并显示错误信息
        fprintf('写入文件 %s 时出错: %s\n', fileName, ME.message);
    end
  
    idxS = find( mask );
    idxN = find( mask2 );
    idxI = find( rotated_mask_dip2 );
    
    % save 
    filename = sprintf('%s/object_def.mat',out_dir);
    % %eval( sprintf('save %s f r dx epsilon_r sigma_e rho idxS idxI idxN mask mask2 dx dy dz pvol',filename) );
    variables_to_save = {'idxS', 'idxI', 'idxN'};
    save(filename, variables_to_save{:});
    
end
