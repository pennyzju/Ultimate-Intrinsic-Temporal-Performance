function rotate_matrix_processing(out_dir, flag, angle_deg)

    % 加载B1m
    load( fullfile(out_dir, 'BASIS_B1m.mat'));
    load( fullfile(out_dir, 'object_def.mat'));

    % 获取数据维度
    m = 95;
    n = 112;
    p = 117;

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
    B1test = zeros(95,112,117);
    B1test(idxS) = B1m(:,1);
    % 旋转矩阵
    rotated_B1test = rotate_matrix(B1test, point1, point2, angle_deg);
    processed_idxS = find(rotated_B1test);

    % save 
    filename = sprintf('%s/object_def.mat',out_dir);
    variables_to_save = {'idxS', 'idxI', 'idxN','processed_idxS'};
    save(filename, variables_to_save{:});
    
end
