function [b1_map, loss_map, usnr_map] = compute_usnr_rf_sta(B1m, rf_all, idxS, idxS0, mask, LOSS)
    % 全 CPU 版本，支持 parfor 并行
    %B1m = , rf_all, idxS, idxS0, mask = usnr_mask_rot;
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
    idx_final = ind_mask(ismember(ind_mask, is_common));  % 输出 map 的索引
    num_vox = length(idx_final);

    % 预分配
    b1_val   = zeros(num_vox, 1);
    loss_val = zeros(num_vox, 1);
    usnr_val = zeros(num_vox, 1);

%     fprintf('ind_mask size: %s\n', mat2str(size(ind_mask)));
% fprintf('idxS size: %s\n', mat2str(size(idxS)));
% fprintf('idxS0 size: %s\n', mat2str(size(idxS0)));

% fprintf('is_common size: %s\n', mat2str(size(is_common)));
% fprintf('loc_mask_in_idxS size: %s\n', mat2str(size(loc_mask_in_idxS)));
% fprintf('loc_idxS_in_mask size: %s\n', mat2str(size(loc_idxS_in_mask)));

% fprintf('loc_mask_in_idxS0 size: %s\n', mat2str(size(loc_mask_in_idxS0)));
% fprintf('loc_idxS0_in_mask size: %s\n', mat2str(size(loc_idxS0_in_mask)));

% fprintf('loc_b1 size: %s\n', mat2str(size(loc_b1)));
% fprintf('loc_rf size: %s\n', mat2str(size(loc_rf)));
% fprintf('idx_final size: %s\n', mat2str(size(idx_final)));
% fprintf('num_vox: %d\n', num_vox);

% fprintf('b1_val size: %s\n', mat2str(size(b1_val)));
% fprintf('loss_val size: %s\n', mat2str(size(loss_val)));
% fprintf('usnr_val size: %s\n', mat2str(size(usnr_val)));


    for i = 1:num_vox

        rf = rf_all(loc_rf(i), :).';         % 列向量
        B1 = B1m(loc_b1(i), :) * rf;
        %loss = real(rf' * LOSS * rf);
        loss= real(rf.' * LOSS * conj(rf));
        usnr = abs(B1) / sqrt(loss);

        b1_val(i)   = B1;
        loss_val(i) = loss;
        usnr_val(i) = usnr;
    end

    b1_map(idx_final)   = b1_val;
    loss_map(idx_final) = loss_val;
    usnr_map(idx_final) = usnr_val;
end
