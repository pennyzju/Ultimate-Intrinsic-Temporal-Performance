function pan_maxit_prepare(public_out_dir, out_dir,fileName, flag, pan_step)
    try
        % 生成文件路径
        filename0 = fullfile(public_out_dir, 'default_material_maps.mat');

        % 检查并创建默认物质映射文件
        if ~exist(filename0, 'file')
            default_material_maps(public_out_dir);
        end

        % 加载默认物质映射数据
        load(filename0);

        % 旋转矩阵
        pan_epsilon_r = shift_gtone_region(epsilon_r, flag, pan_step);
        pan_sigma_e = shift_nonzero_region(sigma_e, flag, pan_step);
        pan_rho = shift_nonzero_region(rho, flag, pan_step);
        pan_mask = double(pan_epsilon_r > 1 & pan_sigma_e > 0 & pan_rho > 0);
        % shave off boundary object voxels
        mask2 = ubasis_shave_boundary_voxels(pan_mask);

        % 将数据写入文件
        save(fullfile(out_dir, fileName), 'dx','pan_epsilon_r', 'pan_sigma_e', 'pan_rho', 'pan_mask','r','mask2');

    catch ME
        % 捕获并显示错误信息
        fprintf('写入文件 %s 时出错: %s\n', fileName, ME.message);
    end

    idxS = find( pan_mask );
    idxN = find( mask2 );
    idxI = find( mask_dip2 );
    

    % save 
    filename = sprintf('%s/object_def.mat',out_dir);
    % %eval( sprintf('save %s f r dx epsilon_r sigma_e rho idxS idxI idxN mask mask2 dx dy dz pvol',filename) );
    variables_to_save = {'idxS', 'idxI', 'idxN'};
    save(filename, variables_to_save{:});
    
end

