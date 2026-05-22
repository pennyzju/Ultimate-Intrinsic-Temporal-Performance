function [sos_map, b1_map, loss_map] = compute_sos_snr_fast(B1m, LOSS, mask, idxS)
    % 输入：
    %   B1m: [N x 2500]，N个已知点的 B1
    %   LOSS: [2500 x 2500]，噪声协方差矩阵
    %   mask: 3D logical matrix，[nx ny nz]
    %   idxS: B1m 每一行对应的线性索引，范围为 1 : (nx*ny*nz)
    % 输出：
    %   sos_map, b1_map, loss_map: 所有值仅在 mask 区域内有值

    [nx, ny, nz] = size(mask);
    N = size(B1m, 1);  % B1m 点数量
    nchannels = size(B1m, 2);

    % ------ Step 1: 建立 idx 映射 -------
    % 找出 mask 中哪些点也存在于 idxS 中
    mask_inds = find(mask);                  % 所有 mask 区域索引
    [common_idx, ia, ib] = intersect(mask_inds, idxS);  % 找到共有点

    % common_idx 是线性索引（在mask中）
    % ia 是 mask_inds 的子索引，ib 是 idxS 的子索引
    if isempty(common_idx)
        warning('mask 区域与 B1m 提供的位置没有交集！');
        sos_map = zeros(nx, ny, nz); b1_map = sos_map; loss_map = sos_map;
        return;
    end

    % ------ Step 2: 提取对应 B1 向量 -------
    B1_selected = B1m(ib, :);  % [M x 2500]

    % ------ Step 3: 计算 SNR -------
    signalmag = real(sum(abs(B1_selected).^2, 2));  % [M x 1]

    B1t = B1_selected.';                    % [2500 x M]
    QB = LOSS * B1t;                        % [2500 x M]
    noisepower = real(sum(conj(B1t) .* QB, 1)).';   % [M x 1]

    sos_snr = signalmag ./ sqrt(noisepower);  % [M x 1]

    % ------ Step 4: 回填结果 -------
    sos_map = zeros(nx, ny, nz);
    sos_map(common_idx) = sos_snr;

    if nargout > 1
        b1_map = zeros(nx, ny, nz);
        b1_map(common_idx) = signalmag;
    end
    if nargout > 2
        loss_map = zeros(nx, ny, nz);
        loss_map(common_idx) = noisepower;
    end
end
