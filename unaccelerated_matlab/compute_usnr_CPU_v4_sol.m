function [usnr_map, b1_map, loss_map, rf_rate, g_factor, rf_all] = compute_usnr_CPU_v4(B1m, LOSS, mask, idxS)

    b1target = 1;
    [nx, ny, nz] = size(mask);
    tmp = zeros(size(mask));
    nfields = size(B1m,2);
    nchannels = nfields;
    b1maps = zeros(nx,ny,nz,nchannels);
    
    % 构建 B1 图
    for i = 1:nfields
        if mod(i,100) == 0
            fprintf('[%d/%d]\n',i,nfields);
        end
        tmp(idxS) = B1m(:,i);
        b1maps(:,:,:,i) = tmp;
    end
    
    q = LOSS;
    
    ind = find(mask > 0);
    [indx, indy, indz] = ind2sub([nx ny nz], ind);
    numInd = numel(ind);
    
    % 预分配结果数组
    tmp1 = zeros(numInd, 1);  % usnr
    tmp2 = zeros(numInd, 1);  % b1
    tmp3 = zeros(numInd, 1);  % loss
    tmp4 = zeros(numInd, 1);  % rf_rate
    tmp5 = zeros(numInd, 1);  % g-factor
    tmp6 = zeros(numInd, nchannels);  % rf_all
    
    startTime = tic;
    parfor_progress(numInd);  % 初始化进度条
    % 前置步骤：构造 b1vectors (numInd × nchannels)，每个点的 B1 向量
    b1vectors = zeros(numInd, nchannels);
    for k = 1:numInd
        b1vectors(k, :) = reshape(b1maps(indx(k), indy(k), indz(k), :), 1, nchannels);
    end
    
    parfor i = 1:numInd
        parfor_progress;
    
        % xpix = indx(i);
        % ypix = indy(i);
        % zpix = indz(i);
    
        % A = compute_system_matrix(b1maps, xpix, ypix, zpix);
        A = b1vectors(i, :);
        A2 = [q, A.'; conj(A), 0];
        b2 = [zeros(nchannels,1); b1target];
    
        %sol = pinv(A2) * b2;
        sol = lsqminnorm(A2, b2);
        rf = conj(sol(1:nchannels));
    
        tmp6(i,:) = rf.';  % 每个点保存一行 rf
    
        tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf;
        loss_ = real(rf.' * q * conj(rf));
    
        tmp1(i) = abs(b1) / sqrt(loss_);
        tmp2(i) = b1;
        tmp3(i) = loss_;
        tmp5(i) = norm(inv(conj(A)*(q \ conj(A'))));
        
        if mod(i, 100) == 0 || i == numInd
            elapsedTime = toc(startTime);
            fprintf('Processed %d of %d (%.2f%%) - Elapsed time: %.2f seconds\n', ...
                i, numInd, (i / numInd) * 100, elapsedTime);
        end
    end
    
    parfor_progress(0);  % 清理进度条
    
    % 构造输出映射
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
        g_factor = zeros(nx, ny, nz);
        g_factor(ind) = tmp5;
    end
    
    if nargout > 5
        rf_all = tmp6;  % 每行是对应 voxel 的 rf 向量
    end
    
    end
    