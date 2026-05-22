function build_new_idxs(B0)
    % 假设 mask 大小为 (95,112,117)
    maskall = zeros(95, 112, 117);  % ✅ 初始化一次

    % 文件夹基础路径
    base_dir = B0_folderpath(B0);
    out_dir = fullfile(base_dir,'left_0');
    flags = {'front','left','rot'};

    for f = 1:length(flags)
        flag = flags{f};

        fprintf('\n正在处理 flag: %s\n', flag);

        for i = -5:1:5
            foldername = fullfile(base_dir, sprintf('%s_%d', flag, i));
            filename = fullfile(foldername, sprintf('%s_%d_material_maps.mat', flag, i));

            if exist(filename, 'file')
                fprintf('  正在处理文件: %s\n', filename);
                data = load(filename);

                if isfield(data, 'rotated_mask')
                    % 🚫 不再判断是否 empty，直接累加
                    maskall = maskall + data.rotated_mask;
                    %fprintf('  当前maskall非零点个数: %d\n', nnz(maskall));
                else
                    warning('文件 %s 中不包含 rotated_mask', filename);
                end
            else
                warning('未找到文件: %s', filename);
            end
        end
    end

    objdef_file = fullfile(base_dir, 'left_0', 'object_def.mat');
    if exist(objdef_file, 'file')
        load(objdef_file, 'idxI', 'idxN');
    else
        error('未找到 object_def.mat 文件: %s', objdef_file);
    end

    idxS = find(maskall > 0);

    filenamenew = fullfile(out_dir, 'object_def_bigmask.mat');
    save(filenamenew, 'idxS', 'idxI', 'idxN', 'maskall');

    fprintf('\n✅ 新的 object_def_bigmask.mat 已保存到: %s\n', filenamenew);
end
