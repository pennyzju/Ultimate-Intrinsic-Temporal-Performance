function [usnr_map, b1_map, loss_map, rf_rate, rf_all] = compute_usnr_gpu_batch(B1m, LOSS, mask, idxS, batchSize)

    if nargin < 5
        batchSize = 512; % 默认批处理大小
    end

    b1target = 1e-3;
    [nx, ny, nz] = size(mask);
    tmp = zeros(size(mask));
    nfields = size(B1m,2);
    nchannels = nfields;

    b1maps = zeros(nx, ny, nz, nchannels);
    for i = 1:nfields
        tmp(idxS) = B1m(:,i);    
        b1maps(:,:,:,i) = tmp;
    end
    q = LOSS;
    ind = find(mask > 0);
    [indx, indy, indz] = ind2sub([nx ny nz], ind);
    numInd = numel(ind);

    tmp1 = zeros(numInd, 1); tmp2 = tmp1; tmp3 = tmp1; tmp4 = tmp1;
    rf_all = zeros(numInd, nchannels);

    startTime = tic;
    for batchStart = 1:batchSize:numInd
        batchEnd = min(batchStart + batchSize - 1, numInd);
        bIdx = batchStart:batchEnd;
        nb = length(bIdx);

        tmp1b = gpuArray.zeros(nb,1);
        tmp2b = gpuArray.zeros(nb,1);
        tmp3b = gpuArray.zeros(nb,1);
        tmp4b = gpuArray.zeros(nb,1);
        rfb = gpuArray.zeros(nb,nchannels);

        q_gpu = gpuArray(q);
        b1target_gpu = gpuArray(b1target);

        for j = 1:nb
            i = bIdx(j);
            xpix = indx(i); ypix = indy(i); zpix = indz(i);
            A = gpuArray(compute_system_matrix(b1maps, xpix, ypix, zpix));
            A2 = [q_gpu, A.'; conj(A), 0];
            b2 = [zeros(nchannels, 1, 'gpuArray'); b1target_gpu];
            sol = pinv(A2) * b2;
            rf = conj(sol(1:nchannels));
            rfb(j,:) = rf.';

            tmp4b(j) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
            b1 = A * rf;
            loss_ = real(rf.' * q_gpu * conj(rf));
            tmp1b(j) = abs(b1) / sqrt(loss_);
            tmp2b(j) = b1;
            tmp3b(j) = loss_;
        end

        % gather 回 CPU
        tmp1(bIdx) = gather(tmp1b);
        tmp2(bIdx) = gather(tmp2b);
        tmp3(bIdx) = gather(tmp3b);
        tmp4(bIdx) = gather(tmp4b);
        rf_all(bIdx,:) = gather(rfb);

        elapsedTime = toc(startTime);
        fprintf('Batch %d-%d processed (%.2f%%) - Elapsed: %.2fs\n', ...
            batchStart, batchEnd, batchEnd/numInd*100, elapsedTime);
    end

    usnr_map = zeros(nx, ny, nz); usnr_map(ind) = tmp1;
    b1_map   = zeros(nx, ny, nz); b1_map(ind)   = tmp2;
    loss_map = zeros(nx, ny, nz); loss_map(ind) = tmp3;
    rf_rate  = zeros(nx, ny, nz); rf_rate(ind)  = tmp4;
end
