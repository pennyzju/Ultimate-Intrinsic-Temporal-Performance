function usnr_mask = create_usnr_mask(flag)
    % 创建用于 USNR 显示的单层掩码
    % 输入:
    %   flag: 'front', 'left', 'rot', 或 'X', 'Y', 'Z'
    % 输出:
    %   usnr_mask: 大小为 [95, 112, 117] 的掩码数组

    % 预设体积大小
    m = 95;
    n = 112;
    p = 117;

    % 初始化掩码
    usnr_mask = zeros(m, n, p);

    % 根据 flag 选择切片方向
    switch flag
        case {'front', 'X'}
            usnr_mask(:, round(n/2), :) = 1;
        case {'left', 'Z'}
            usnr_mask(round(m/2), :, :) = 1;
        case {'rot', 'Y'}
            usnr_mask(:, :, round(p*0.74)) = 1;  % 原来是87/117 ≈ 0.74
        otherwise
            error('未知的 flag 值，请使用 "front"、"left"、"rot" 或 "X"、"Y"、"Z"');
    end
end
