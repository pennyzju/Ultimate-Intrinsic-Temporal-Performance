function ubasis_solve_new(public_out_dir,out_dir,fileName)

    load(fullfile(public_out_dir, 'ubasis_options.mat'));
    load(fullfile(out_dir, fileName));
    
    filename = sprintf('%s/object_def.mat',out_dir);
    load(filename);
    
    log_filename = sprintf('%s/LOG_solve.txt',out_dir);
    fid = fopen(log_filename,'w');
    
    % -------------------------------------------------------------------------
    %                  EM constants
    % -------------------------------------------------------------------------
    gamma = 42.576e6;
    
    TESLA = f/gamma;
    mu = 4*pi*1e-7;
    co = 299792458;
    eo = 1/co^2/mu;
    
    omega = 2 * pi * f;
    lambda  = co/f;
    ko = 2*pi/lambda;
    eta =  3.767303134617706e+002; % Free-space impedance

    e_r = rotated_epsilon_r - 1j*rotated_sigma_e/(eo*omega);%注意原始矩阵的大小
    
    % -------------------------------------------------------------------------
    %            DATA
    % -------------------------------------------------------------------------
    
    
    [L,M,N] = size(e_r);
    nD = L*M*N;
    nS = length(idxS);
    nI = length(idxI);
    dV = dx^3;
    [~, maxeps] = max(real(e_r(:)));
    
    idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid
    idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid
    
    
    fprintf(fid, '\n\n ----------------------------------------------------------');
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n');
    fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
    fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
    fprintf(fid, '\n     # DOFS:                %d', 3*L*M*N);
    fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
    fprintf(fid, '\n     # DOFS in Coil Shell:  %d', 3*length(idxI));
    fprintf(fid, '\n     Operating Frequency:   %.2f MHz (%.1f T)',f/1e6,f/(gamma));
    fprintf(fid, '\n     Maximum Epsilon:       %.2f %.2f j', real(e_r(maxeps)), imag(e_r(maxeps)));
    fprintf(fid, '\n');
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n ----------------------------------------------------------\n\n ');
    type(log_filename);
    
    
    % -------------------------------------------------------------------------
    %                  Prepare the cases to Solve
    % -------------------------------------------------------------------------
    
    % define subdomain of interest
    xd = r(:,:,:,1);
    xmin = min(xd(idxS))-dx; xmax = max(xd(idxS))+dx;
    yd = r(:,:,:,2);
    ymin = min(yd(idxS))-dx; ymax = max(yd(idxS))+dx;
    zd = r(:,:,:,3);
    zmin = min(zd(idxS))-dx; zmax = max(zd(idxS))+dx;
    
    % cut grid to subdomain dimensions
    idxX = find( (r(:,1,1,1) >= xmin) & (r(:,1,1,1) <= xmax));
    idxY = find( (r(1,:,1,2) >= ymin) & (r(1,:,1,2) <= ymax));
    idxZ = find( (r(1,1,:,3) >= zmin) & (r(1,1,:,3) <= zmax));
    
    
    % -------------------------------------------------------------------------
    %                  Create new grid
    % -------------------------------------------------------------------------
    
    % truncate e_r to new reduced grid
    e_rnew = e_r(idxX,idxY,idxZ);
    rho_new = rotated_rho(idxX,idxY,idxZ);
    
    % get dimensions
    [Lnew,Mnew,Nnew] = size(e_rnew);
    
    % -------------------------------------------------------------------------
    %                  Determine good circulant size and compute if needed
    % -------------------------------------------------------------------------
    
    
    % get the required circulant for the EF operator
    [fG] = retrieve_fG(Lnew,Mnew,Nnew,dx,f);
    
    
    % -------------------------------------------------------------------------
    %                 Define the solution parameters
    % -------------------------------------------------------------------------
    
    % prepare domain
    idxE = find( e_rnew>1 & rho_new> 0 );
    idxE3 = [idxE; Lnew*Mnew*Nnew+idxE; 2*Lnew*Mnew*Nnew+idxE]; % the vector of non-air positions in 3D grid
    
    infocir = whos('fG');
    
    if 1  % set to 0 if memory management problem during the MARIE solve
        
        dev = gpuDevice;
        reset(dev);                
        
        % define the main function
        fprintf(fid,'\n Estimated MVP Peak Memory %.3f MB.\n Applying GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fprintf('\n Estimated MVP Peak Memory %.3f MB.\n Applying GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        %type(log_filename);        
        
        % fA   = @(J,transp_flag)mv_op_fft_Nop_2_gpu_bastien(J, fG, e_rnew, dV, transp_flag, idxE3);
        
        FG_gpu = gpuArray(fG);
        EPS_gpu = gpuArray(e_rnew);
        fA   = @(J,transp_flag)mv_op_fft_Nop_2_gpu(J, FG_gpu, EPS_gpu, dV, transp_flag, idxE3);
    
    else
        fprintf(fid,'\n Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fprintf('\n Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        %type(log_filename);
        fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_rnew, dV, sTrans, idxE3);
            
    end
        
    
    % -------------------------------------------------------------------------
    % Obtain the Scattered field due to Einc
    % -------------------------------------------------------------------------
    
    % set the multiplier to convert fields to suitable input
    tau = 1j*omega*eo*(e_rnew(:)-1)./(e_rnew(:));
    mult = [tau; tau; tau];
    mult = dV*mult(idxE3); % chose Scatterer positions
    
    
    % Load the basis
    Out_file = sprintf('%s/BASIS_Einc.mat', out_dir);
    load(Out_file);
    Jsol = Ebasis; clear Ebasis;
    ncol = size(Jsol,2);
    
    fprintf(fid, '\n  Solving the system for %d inputs \n', ncol);
    fprintf('\n  Solving the system for %d inputs \n', ncol);
    %type(log_filename);
    t1 = clock;
    
    INDEX_SAVE_PRE = 1;
    
    for ii=INDEX_SAVE_PRE:ncol
        
        % get the incident field excitation and transform into currents
        JExc = mult.*Jsol(:,ii);
        
        % solve to obtain the total currents
        [jsol,flag,relres,iter,resvec] = bicgstab(@(J)fA(J,'notransp'), JExc, tol, maxit);
        
        if flag~= 0        
            warning( sprintf('\t\tbicgstab did not converge: FLAG = %3.0f   ITER = %3.1f\n',flag,iter) );
        end
        
       Jsol(:,ii) = jsol; % store the solution currents
        
        if 1  % ((ii/100) == floor(ii/100))
            t2 = clock;
            fprintf(fid, '\n %d Solves: Elapsed time = %f   flag = %3.0f   iter = %3.1f   relres = %e', ii, etime(t2,t1), flag, iter, relres);
            fprintf('\n %d Solves: Elapsed time = %f   flag = %3.0f   iter = %3.1f   relres = %e', ii, etime(t2,t1), flag, iter, relres);
            %type(log_filename);
        end
        
        % SAVE JSOL EVERY 100 FIELDS IN CASE OF SIMULATION CRASH
        if mod(ii,100) == 0 || ii==ncol
            
            TEMP = Jsol( :,INDEX_SAVE_PRE :ii );
            eval( sprintf('save -v7.3 Jsol_TEMP_%dTO%d.mat TEMP;',INDEX_SAVE_PRE,ii) );
            
            INDEX_SAVE_PRE = ii+1;        
        end
            
        
    end
    
    clear jsol; clear JExc;
            
    t2 = clock;
    fprintf(fid, '\n\n Equivalent current generated for %d vectors. Elapsed time %f \n', ncol, etime(t2,t1));
    fprintf(fid, '\n\n Equivalent current generated for %d vectors. Elapsed time %f \n', ncol, etime(t2,t1));
    %type(log_filename);
        
    Out_file = sprintf('%s/BASIS_Jsol.mat', out_dir);
    save(Out_file , 'Jsol', '-v7.3');
    
    fclose(fid);
    
    
    
    
    
    