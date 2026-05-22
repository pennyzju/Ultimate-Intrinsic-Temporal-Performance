function process_usnr_main(B0, flag, angle_deg)
    % 封装 UISNR 主流程
    % 输入：
    %   B0 - 磁场强度（支持 7、3、1.5）
    %   flag - 'front'、'left'、'rot' 之一
    %   angle_deg - 旋转角度（整数）
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/MARIE_ubasis_v6_noSVD_BIGMEM'));
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/unaccelerated_matlab'));
    addpath(genpath('/data/jiaxinli/projects/20240728_my-repo_test/RUN_DUKE_BASIS_SET'));

    if B0 == 7
        FolderPath = '/data2/jiaxinli/pan/20240801_UISNR_output';
    elseif B0 == 3
        FolderPath = '/data2/jiaxinli/pan/20241101_3T_UISNR_output';
    elseif B0 == 1.5
        FolderPath = '/data2/jiaxinli/pan/20250201_1_5T_UISNR_output';
    else
        error('B0 must be 7, 3, or 1.5');
    end

    folder_name = fullfile(FolderPath,'usnr_check');
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);  % 如果文件夹不存在，则创建它
    end

    subfoldername = sprintf('%s_%d', flag, angle_deg);
    target_path = fullfile(FolderPath, subfoldername);
    if exist(target_path, 'dir')
        cd(target_path);
        load BASIS_B1m.mat;
        load('object_def.mat', 'idxS');
        load BASIS_LOSS_N.mat;
    else
        error('目录不存在：%s', target_path);
    end

    [point1, point2] = get_rotation_axis(flag);
    %load ../public/mask_3plain;
    mask = create_usnr_mask(flag);
    %mask = rotate_matrix(mask, point1, point2, angle_deg);
    save(fullfile(folder_name, sprintf('mask_%s.mat', flag)), 'mask');

    flag = [1 2]; 
    %1 - E+M(1:1)
    %2 - E OR M
    %3 - 0.75 E AND 0.25 M
    %3 - 0.25 E AND 0.75 M

    if ismember(1,flag) == 1
        %E+M(1:1)
        seq = [500:1000:2500];
        for i=1:length(seq)
            % if exist(['../SNR/USNR_',num2str(seq(i)),'_full.mat'],'file')==0
            %     if seq(i)<100
            %         load(['BASIS_B1m_100.mat']);
            %         B1m = B1m(:,1:seq(i));
            %     else
            %         load(['BASIS_B1m_',num2str(seq(i)),'.mat']);
            %     end
            %load mask_plain.mask;
            LOSS_seq=LOSS(1:seq(i),1:seq(i));
            B1m_seq = B1m(1:seq(i));
            [nx,ny,nz] = size(mask);
            % nfields = i;
            % nchannels = nfields;
            [usnr_map,b1_map,loss_map]=compute_usnr_unacc_fast_v3_GAO_noParallel(B1m_seq,LOSS_seq,mask,idxS);
            out_snr_file = fullfile(folder_name, sprintf('USNR_%s_%d_full.mat', flag, angle_deg));
            save(out_snr_file, 'usnr_map', '-v7.3');
            out_loss_file = fullfile(folder_name, sprintf('loss_%s_%d_full.mat', flag, angle_deg));
            save(out_snr_file, 'usnr_map', '-v7.3');
            clear usnr_map b1_map loss_map B1m_seq LOSS_seq
        end
    end
end