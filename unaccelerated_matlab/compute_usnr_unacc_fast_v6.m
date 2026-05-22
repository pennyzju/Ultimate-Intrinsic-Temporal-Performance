function [usnr_map, b1_map, loss_map, rf_rate, rf_map] = compute_usnr_unacc_fast_v6(B1m, LOSS, mask, idxS)
    % GPU + 矢量化并行版本，加速 USNR / B1 / Loss / rf_rate 计算
    
    b1target = 1e-3;
    [nx, ny, nz] = size(mask);
    nfields = size(B1m, 2);
    nchannels = nfields;
    
    % 构建 4D b1maps（nx, ny, nz, nchannels）
    tmp = zeros(size(mask));
    b1maps = zeros(nx, ny, nz, nchannels);
    for i = 1:nfields
        tmp(idxS) = B1m(:, i);
        b1maps(:, :, :, i) = tmp;
    end
    
    % 获取 mask 中的体素索引（线性和三维）
    ind = find(mask > 0);
    [indx, indy, indz] = ind2sub([nx, ny, nz], ind);
    numInd = numel(ind);
    
    % 提取所有 mask 中的 A 向量（numInd × nchannels）
    b1maps_reshaped = reshape(b1maps, [], nchannels);
    A_all = b1maps_reshaped(ind, :);  % double[numInd x nchannels]
    
    % 转移到 GPU
    A_all = gpuArray(A_all);
    q_gpu = gpuArray(LOSS);
    b1target_gpu = gpuArray(b1target);
    
    % 快速计算 rf_all（numInd × nchannels）
    AQ = A_all / q_gpu;                             % numInd × nchannels
    denom = sum(AQ .* conj(A_all), 2);              % numInd × 1
    rf_all = conj((AQ') .* (b1target_gpu ./ denom.'));  % nchannels × numInd
    rf_all = rf_all.';                              % numInd × nchannels
    
    % B1 计算（numInd × 1）
    b1 = sum(A_all .* rf_all, 2);
    
    % Loss 计算（numInd × 1）
    loss = real(sum((rf_all * q_gpu) .* conj(rf_all), 2));
    
    % USNR 计算
    usnr = abs(b1) ./ sqrt(loss);
    
    % rf_rate 计算
    odd_idx = 1:2:nchannels;
    even_idx = 2:2:nchannels;
    rf_rate = sum(abs(rf_all(:, odd_idx)), 2) ./ sum(abs(rf_all(:, even_idx)), 2);
    
    % 将结果 gather 回 CPU
    usnr = gather(usnr);
    b1 = gather(b1);
    loss = gather(loss);
    rf_rate = gather(rf_rate);
    rf_all = gather(rf_all);
    
    % 写入输出图
    usnr_map = zeros(nx, ny, nz);
    usnr_map(ind) = usnr;
    
    if nargout > 1
        b1_map = zeros(nx, ny, nz);
        b1_map(ind) = b1;
    end
    if nargout > 2
        loss_map = zeros(nx, ny, nz);
        loss_map(ind) = loss;
    end
    if nargout > 3
        rf_rate_map = zeros(nx, ny, nz);
        rf_rate_map(ind) = rf_rate;
        rf_rate = rf_rate_map;
    end
    if nargout > 4
        %rf_map = zeros(nx, ny, nz);
        rf_map = rf_all;
    end
    end
    