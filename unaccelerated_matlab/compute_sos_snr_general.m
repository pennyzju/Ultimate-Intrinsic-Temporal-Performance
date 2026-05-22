function [sos_snr_map, b1_map, loss_map, rf_rate] = compute_sos_snr_general(B1m, q, mask, idxS, abs_rf, abs_q)
    % 通用计算 SNR 的函数
    % 输入:
    %   B1m: [Npoints x Nchannels] B1矩阵
    %   q: 噪声协方差矩阵 [Nchannels x Nchannels]
    %   mask: 3D掩码数组
    %   idxS: mask对应的线性索引，用于填充3D矩阵
    %   opts: 结构体，可选参数：
    %       opts.abs_rf (bool): 是否对 rf 取绝对值，默认 false
    %       opts.abs_q (bool): 是否对 q 取绝对值，默认 false
    %
    % 输出:
    %   sos_snr_map, b1_map, loss_map, rf_rate 均为3D矩阵
    
    arguments
        B1m (:,:) double
        q (:,:) double
        mask (:,:,:) double
        idxS (:,1) double
        abs_rf (1,1) logical = false
        abs_q (1,1) logical = false
    end

    nfields = size(B1m, 2);
    nchannels = nfields;
    [nx, ny, nz] = size(mask);
    b1maps = zeros(nx, ny, nz, nchannels);
    tmp = zeros(size(mask));

    % 填充b1maps
    for i = 1:nfields
        if mod(i, 100) == 0
            fprintf('[%d/%d]\n', i, nfields);
        end
        tmp(idxS) = B1m(:, i);
        b1maps(:,:,:,i) = tmp;
    end

    ind = find(mask > 0);
    [indx, indy, indz] = ind2sub([nx, ny, nz], ind);

    % 初始化输出矩阵
    b1_map = zeros(nx, ny, nz);
    sos_snr_map = zeros(nx, ny, nz);
    loss_map = zeros(nx, ny, nz);
    rf_rate = zeros(nx, ny, nz);

    numInd = numel(ind);

    % 预分配中间变量
    tmp_snr = zeros(numInd, 1);
    tmp_b1 = zeros(numInd, 1);
    tmp_loss = zeros(numInd, 1);
    tmp_rf_rate = zeros(numInd, 1);

    startTime = tic;

    q_mat = q;
        if abs_q
            q_mat = abs(q_mat);
        end

    for i = 1:numInd
        xpix = indx(i);
        ypix = indy(i);
        zpix = indz(i);

        % 提取该点所有通道的 B1 值
        A = compute_system_matrix(b1maps, xpix, ypix, zpix);

        b1 = norm(A); % ||A|| = sqrt(sum(|A_i|^2))

        rf = conj(A(1:nchannels));
        if abs_rf
            rf = abs(rf);
        end

        rf = rf(:);       
        loss_= (rf.')*q*conj(rf) ;
        %loss_ = real(rf.' * q_mat * conj(rf));

        tmp_snr(i) = abs(b1) / sqrt(loss_);
        tmp_b1(i) = b1;
        tmp_loss(i) = loss_;
        tmp_rf_rate(i) = 0; % 如果需要计算rf_rate，后续可添加逻辑

        if mod(i, 100) == 0 || i == numInd
            elapsedTime = toc(startTime);
            fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f seconds\n', ...
                i, numInd, (i / numInd) * 100, elapsedTime);
        end
    end

    sos_snr_map(ind) = tmp_snr;
    b1_map(ind) = tmp_b1;
    loss_map(ind) = tmp_loss;
    rf_rate(ind) = tmp_rf_rate;
end
