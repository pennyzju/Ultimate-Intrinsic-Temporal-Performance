function [usnr_map, b1_map, loss_map, rf_rate, rf_all] = compute_usnr_unacc_fast_v5(B1m, LOSS, mask, idxS)
    %与V4相比增加了输出rf功能
    b1target = 1e-3;
    [nx, ny, nz] = size(mask);

    tmp = zeros(size(mask));
    nfields = size(B1m,2);  % 通道数
    nchannels = nfields;

    % 构建 4D B1map [nx, ny, nz, nchannels]
    b1maps = zeros(nx, ny, nz, nchannels);
    for i = 1:nfields
        if mod(i,100) == 0
            fprintf('[%d/%d]\n',i,nfields);
        end
        tmp(idxS) = B1m(:,i);    
        b1maps(:,:,:,i) = tmp;
    end
    clear tmp;

    q = LOSS;

    ind = find(mask > 0);
    [indx, indy, indz] = ind2sub([nx ny nz], ind);
    numInd = numel(ind);

    % 初始化 GPU 上的临时变量
    tmp4 = gpuArray.zeros(numInd, 1);  % rf_rate
    tmp1 = gpuArray.zeros(numInd, 1);  % usnr
    tmp2 = gpuArray.zeros(numInd, 1);  % b1
    tmp3 = gpuArray.zeros(numInd, 1);  % loss
    rf_all_gpu = gpuArray.zeros(numInd, nchannels);  % 新增：保存 rf 向量

    q_gpu = gpuArray(q);
    b1target_gpu = gpuArray(b1target);
    indx_gpu = gpuArray(indx);
    indy_gpu = gpuArray(indy);
    indz_gpu = gpuArray(indz);

    startTime = tic;
    for i = 1:numInd
        xpix = indx_gpu(i);
        ypix = indy_gpu(i);
        zpix = indz_gpu(i);

        A = gpuArray(compute_system_matrix(b1maps, xpix, ypix, zpix));
        % fprintf('size(q_gpu) = [%d, %d]\n', size(q_gpu));
        % fprintf('size(A.'') = [%d, %d]\n', size(A.'));
        % fprintf('size(conj(A)) = [%d, %d]\n', size(conj(A)));
        A2 = [q_gpu, A.'; conj(A), 0];
        b2 = [zeros(nchannels, 1, 'gpuArray'); b1target_gpu];

        sol = pinv(A2) * b2;
        rf = conj(sol(1:nchannels));

        rf_all_gpu(i, :) = rf.';  % 保存 rf（行向量）

        tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf;
        loss_ = real(rf.' * q_gpu * conj(rf));

        tmp1(i) = abs(b1) / sqrt(loss_);
        tmp2(i) = b1;
        tmp3(i) = loss_;

        if mod(i, 100) == 0 || i == numInd
            elapsedTime = toc(startTime);
            fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f seconds\n', ...
                i, numInd, (i / numInd) * 100, elapsedTime);
        end
    end

    % 将结果从 GPU 拷回 CPU
    tmp1 = gather(tmp1);
    tmp2 = gather(tmp2);
    tmp3 = gather(tmp3);
    tmp4 = gather(tmp4);
    rf_all = gather(rf_all_gpu);

    usnr_map = zeros(nx, ny, nz);
    usnr_map(ind) = tmp1;

    if nargout > 1
        b1_map = zeros(nx, ny, nz);
        b1_map(ind) = tmp2;
    end

    if nargout > 2
        loss_map = zeros(nx, ny, nz);
        loss_map(ind) = tmp3;
    end

    if nargout > 3
        rf_rate = zeros(nx, ny, nz);
        rf_rate(ind) = tmp4;
    end

    if nargout > 4
        % rf_all 已在循环中构造并返回
    end
end
