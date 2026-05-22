function ubasis_comp_fields_new_mag(public_out_dir,out_dir,fileName,subfolder)
    ubasis_options = load(fullfile(public_out_dir, 'ubasis_options.mat'));
    
    object_def = load(fullfile(out_dir, 'object_def.mat'));
    load(fullfile(out_dir, fileName));%导入材料数据
    
    % load basis solution
    BASIS_Jsol =load(fullfile(out_dir, 'BASIS_Jsol.mat')); 
    Nexc = size(BASIS_Jsol.Jsol,2);% store both the E and H dipoles (2*Nexc)
    
    log_filename = sprintf('%s/LOG_comp_fields.txt',out_dir);
    fid = fopen(log_filename,'w');
    
    %%  DATA
    
    [L,M,N] = size(rotated_epsilon_r);
    dV = dx^3;
    pvol = dx^3;
    nD = L*M*N;
    nS = length(object_def.idxS);
    nN = length(object_def.idxN);

    idxS3 = [object_def.idxS; nD+object_def.idxS; 2*nD+object_def.idxS]; % the vector of non-air positions in 3D grid
    idxN3 = [object_def.idxN; nD+object_def.idxN; 2*nD+object_def.idxN]; % the vector of noise related positions in 3D grid
    
    %%  type(log_filename);
    fprintf(fid, '\n\n ----------------------------------------------------------');
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n');
    fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
    fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
    fprintf(fid, '\n     # DOFS:                %d', 3*L*M*N);
    fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*nS);
    fprintf(fid, '\n     # DOFS in Coil Shell:  %d', 3*nN);
    fprintf(fid, '\n');
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n ----------------------------------------------------------\n\n ');
        
    type(log_filename);
    %%  EM constants
    % -------------------------------------------------------------------------
    
    %gamma = 42.576e6;
    
    %Tesla = ubasis_options.f/gamma;
    mu = 4*pi*1e-7;
    co = 299792458;
    eo = 1/co^2/mu;
    %
    omega = 2 * pi * ubasis_options.f;
    %lambda  = co/ubasis_options.f;
    %ko = 2*pi/lambda;
    %omega_mu = omega*mu;
    %eta =  3.767303134617706e+002; % Free-space impedance
    
    % e_r = epsilon_r - 1j*sigma_e/(eo*omega);
    
    %%  get the Circulants for the radiated fields
    ttot = tic;
    
    % operators with LMN original size
    % get the required circulant for the EF operator
    [fN] = retrieve_fG(L,M,N,dx,ubasis_options.f);
    
    timecirc = toc(ttot);
    fprintf(fid,'\n\n Circulants generated, time %f', timecirc);
    fprintf('\n\n Circulants generated, time %f', timecirc);
    %type(log_filename);
    
    % -------------------------------------------------------------------------
    % Apply the direct operator on the FN case
    % -------------------------------------------------------------------------
    
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n N Operator starting with %d excitations \n', Nexc);
    fprintf('\n\n N Operator starting with %d excitations \n', Nexc);
    %type(log_filename);
    
    % GPU use
    dev = gpuDevice(1);
    reset(dev);
    
    % load basis
    Out_file = sprintf('%s/BASIS_Einc.mat', out_dir);
    BASIS_Einc = load(Out_file);
    
    % allocate space
    ve = zeros(L,M,N,3);
    vj = zeros(L,M,N,3);
    Nexc = size(BASIS_Einc.Ebasis,2);
    
    Ebasis = zeros(3*nS,Nexc);
    Enoise = zeros(3*nN,Nexc);
    
    t1 = tic;
    t2 = tic;
    for ii = 1:Nexc
        
        ve(idxS3) = BASIS_Einc.Ebasis(:,ii); % transfer to global vars
        vj(idxS3) = BASIS_Jsol.Jsol(:,ii); % transfer to global vars
        
        % compute total field radiated
        [Ein] = E_field_Nop_compGPU_bastien(vj, fN, dV, omega, eo, ve);
        
        % transfer to local vars
        Ebasis(:,ii) = Ein(idxS3);
        Enoise(:,ii) = Ein(idxN3);
        
        if ((ii/100) == floor(ii/100))
            fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
            fprintf('\n %d Fields. Elapsed time %g', ii, toc(t2));
            %type(log_filename);
            t2 = tic;
        end
        
    end
    
    
    % FOR HYPERTHERMIA ONLY!!!
    save -v7.3 BASIS_Etot.mat Ebasis;
    
    % clear FNgpu;
    reset(dev);
    clear Ein; clear ve; clear vj;
    clear fN;
    
    timeinc = toc(t1);
    fprintf(fid,'\n\n Total fields generated generated, time %f', timeinc);
    fprintf('\n\n Total fields generated generated, time %f', timeinc);
    
    %% compute the loss and gsar matrices
    t1 = tic;
    
    Ebasis = reshape(Ebasis,nS,3,Nexc); % reshape in x y and z components
    Enoise = reshape(Enoise,nN,3,Nexc); % reshape in x y and z components
    
    Sig = rotated_sigma_e(object_def.idxS);
    Rho = rotated_rho(object_def.idxS);
    
    % 计算每一点电场模值平方（nS × Nexc）
    Emag_sq = squeeze(sum(abs(Ebasis).^2, 2));  % [nS × Nexc]

    LOSS = (Emag_sq.' * (Sig * pvol));  % Nexc × 1
    
    % Emag_sq: [nS × Nexc]，每个激发在每个体素上的 |E|²
    Emag_sq = squeeze(sum(abs(Ebasis).^2, 2));  % nS × Nexc

    % 加权矩阵：W_diag 是 nS × nS 的稀疏对角阵（权重 = σ * pvol）
    W_diag = spdiags(Sig * pvol, 0, nS, nS);  % 加入体素体积

    % 最终 LOSS 为 Nexc × Nexc，对应每两个 basis 之间的功率耦合
    LOSS = Emag_sq' * W_diag * Emag_sq;  % (Nexc × nS) * (nS × nS) * (nS × Nexc) = [Nexc × Nexc]

    % 确保结果是实数（由于数值误差，可能出现微小虚部）
    LOSS = real(LOSS);  % 单位：瓦特（W）

    % RHODIAG = spdiags(Sig./Rho,0,nS,nS);
    % GSAR = (squeeze(Ebasis(:,1,:)))'*(RHODIAG*(squeeze(Ebasis(:,1,:)))) + (squeeze(Ebasis(:,2,:)))'*(RHODIAG*(squeeze(Ebasis(:,2,:)))) + (squeeze(Ebasis(:,3,:)))'*(RHODIAG*(squeeze(Ebasis(:,3,:))));
    
    % npix = size(object_def.idxS,1);
    % GSAR = GSAR.' / (2*npix);  % SAR is in W/kg averaged over the whole volume (hence the factor npix)
    %计算 GSAR: [Nexc × 1]
    GSAR = Emag_sq.' * (Sig ./ Rho);              % (Nexc × nS) * (nS × 1) → [Nexc × 1]
    GSAR = GSAR / (2 * nS);
    
    clear Ebasis;
    
    t2 = toc(t1);
    fprintf(fid, '\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
    fprintf('\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
    
    % type(log_filename);
    output_path = fullfile(out_dir, subfolder);
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    Xout_file = fullfile(output_path, 'BASIS_LOSS.mat');
    save(Xout_file, 'LOSS', '-v7.3');

    Xout_file = fullfile(output_path, 'BASIS_GSAR.mat');
    save(Xout_file, 'GSAR', '-v7.3');
    
    t1 = tic;
    
    SigN = rotated_sigma_e(object_def.idxN);    % [nN × 1]
    RhoN = rotated_rho(object_def.idxN);        % [nN × 1]

    w_loss_N = SigN * pvol;                 % [nN × 1]
    w_gsar_N = (SigN ./ RhoN);              % [nN × 1]

    % 电场 reshape: [nN × 3 × Nexc]
    Enoise = reshape(Enoise, nN, 3, Nexc);

    % 电场模值平方: [nN × Nexc]
    Emag_sq_N = sum(abs(Enoise).^2, 2);     % [nN × 1 × Nexc]
    Emag_sq_N = squeeze(Emag_sq_N);         % [nN × Nexc]

    % 加权后的电场幅值：每个像素乘 sqrt(σ·dv)，类似 sqrt 权重
    Ew = Emag_sq_N .* sqrt(w_loss_N);          % [nN × Nexc]

    % 计算 LOSS 矩阵（Nexc × Nexc）——近似耦合矩阵，仅幅值近似
    LOSS = Ew' * Ew;                        % [2500 × 2500]，单位 W

    % 计算 GSAR_N: [Nexc × 1]
    GSAR = Emag_sq_N.' * w_gsar_N;
    GSAR = GSAR / (2 * nN);

    
    t2 = toc(t1);
    fprintf(fid, '\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
    fprintf('\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
    % type(log_filename);
    
    Xout_file = fullfile(output_path, 'BASIS_LOSS_N.mat');
    save(Xout_file, 'LOSS', '-v7.3');

    Xout_file = fullfile(output_path, 'BASIS_GSAR_N.mat');
    save(Xout_file, 'GSAR', '-v7.3');
    
    % -------------------------------------------------------------------------
    %            get the Circulants for the radiated fields
    % -------------------------------------------------------------------------
    
    ttot = tic;
    
    % operators with LMN original size
    
    % get the required circulant for the MF operator
    [fK] = retrieve_fK(L,M,N,dx,ubasis_options.f);
    
    timecirc = toc(ttot);
    fprintf(fid,'\n\n Circulants generated, time %f', timecirc);
    fprintf('\n\n Circulants generated, time %f', timecirc);
    % type(log_filename);
    
    % -------------------------------------------------------------------------
    % Apply the direct operator on the FK case
    % -------------------------------------------------------------------------
    
    
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n K Operator starting with %d excitations \n', Nexc);
    fprintf('\n\n K Operator starting with %d excitations \n', Nexc);
    %type(log_filename);
    
    % GPU use
    dev = gpuDevice(1);
    reset(dev);
    % send data to GPU mem
    % FKgpu = gpuArray(fK);
    
    % load basis
    Out_file = sprintf('%s/BASIS_Hinc.mat', out_dir);
    BASIS_Hinc = load(Out_file);
    
    % allocate space
    vm = zeros(L,M,N,3);
    vj = zeros(L,M,N,3);
    Nexc = size(BASIS_Hinc.Hbasis,2);
    
    t1 = tic;
    t2 = tic;
    Hbasis = zeros(3*nS,Nexc);
    
    for ii = 1:Nexc
        
        vm(idxS3) = BASIS_Hinc.Hbasis(:,ii); % transfer to global vars
        vj(idxS3) = BASIS_Jsol.Jsol(:,ii); % transfer to global vars
        
        % compute total field radiated
        % incident field is zero (0*ve)
        [Hin] = H_field_Kop_compGPU_bastien(vj, fK, dV, vm);
        
        % transfer to local vars
        Hbasis(:,ii) = mu*Hin(idxS3);
        
        if ((ii/500) == floor(ii/500))
            fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
            fprintf('\n %d Fields. Elapsed time %g', ii, toc(t2));
            % type(log_filename);
            t2 = tic;
        end
        
    end
    
    % clear FKgpu; 
    reset(dev);
    clear vj; clear ve; clear Hin;
    clear fK;
    
    timeinc = toc(t1);
    fprintf(fid,'\n\n Total fields generated generated, time %f', timeinc);
    fprintf('\n\n Total fields generated generated, time %f', timeinc);
    % type(log_filename);
    
    % -------------------------------------------------------------------------
    %            compute the B1p and B1m matrices
    % -------------------------------------------------------------------------
    
    t1 = tic;
    
    Hbasis = reshape(Hbasis,nS,3,Nexc); % reshape in x y and z components
    
    B1p = ( squeeze(Hbasis(:,1,:)) + 1j*squeeze(Hbasis(:,2,:)) ) / sqrt(2);
    Xout_file = fullfile(output_path, 'BASIS_B1p.mat');
    save(Xout_file, 'B1p', '-v7.3');
    
    B1m = ( squeeze(Hbasis(:,1,:)) - 1j*squeeze(Hbasis(:,2,:)) ) / sqrt(2);
    Xout_file = fullfile(output_path, 'BASIS_B1m.mat');
    save(Xout_file, 'B1m', '-v7.3');
    
    clear Hbasis;
    clear B1p; clear B1m;
    
    t2 = toc(t1);
    fprintf(fid, '\n\n B1p and B1m matrices computed. Elapsed time %f \n', t2);
    fprintf('\n\n B1p and B1m matrices computed. Elapsed time %f \n', t2);
    % type(log_filename);
    
    fclose(fid);
    
    
    
    
    
    
    
    
    
    