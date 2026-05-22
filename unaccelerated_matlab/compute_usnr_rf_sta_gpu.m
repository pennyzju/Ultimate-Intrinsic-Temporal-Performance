function [b1_map, loss_map, usnr_map] = compute_usnr_rf_sta_gpu(B1m, rf_all, idxS, idxS0, mask, LOSS)
    % 使用 GPU 单线程计算，无 parfor，避免 GPU 冲突

    [nx, ny, nz] = size(mask);
    b1_map   = zeros(nx, ny, nz);
    loss_map = zeros(nx, ny, nz);
    usnr_map = zeros(nx, ny, nz);

    % 提取 mask 中 >0 的索引
    ind_mask = find(mask(:) > 0);

    % 找出这些点在 B1 和 rf 中的对应行
    [is_common, loc_mask_in_idxS, loc_idxS_in_mask] = intersect(idxS, ind_mask);
    [~, loc_mask_in_idxS0, loc_idxS0_in_mask] = intersect(idxS0, ind_mask);

    loc_b1 = loc_mask_in_idxS;
    loc_rf = loc_mask_in_idxS0;
    idx_final = ind_mask(ismember(ind_mask, is_common));
    num_vox = length(idx_final);

    % GPU 数据
    B1m_gpu = gpuArray(B1m);
    rf_all_gpu = gpuArray(rf_all);
    LOSS_gpu = gpuArray(LOSS);

    % 初始化
    b1_val   = zeros(num_vox, 1);
    loss_val = zeros(num_vox, 1);
    usnr_val = zeros(num_vox, 1);

    for i = 1:num_vox
        rf = rf_all_gpu(loc_rf(i), :).';     % 列向量
        B1 = B1m_gpu(loc_b1(i), :) * rf;
        loss = real(rf' * LOSS_gpu * rf);
        usnr = abs(B1) / sqrt(loss + eps);

        b1_val(i)   = gather(B1);
        loss_val(i) = gather(loss);
        usnr_val(i) = gather(usnr);
    end

    % 写入三维图
    b1_map(idx_final)   = b1_val;
    loss_map(idx_final) = loss_val;
    usnr_map(idx_final) = usnr_val;
end
