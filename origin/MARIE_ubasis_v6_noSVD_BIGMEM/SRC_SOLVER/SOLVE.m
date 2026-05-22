function [JSOL,FLAG,RELRES,ITER,RESVEC,TIME] = SOLVE(EINC,FG,EPS,R,FREQ,FORMULATION,LOCAL,METHOD,TOL,IT,INNER_IT,OUTER_IT,RITZ,GPU,PRECOND,LOGFILE)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%
%
%   INPUT:  EINC        Incident field
%           FG          Circulant matrix for the given geometry, freq and resolution
%           EPS         3D matrix with EM characteristics of the non-air media
%           R           3D vector with the domain grid
%           FREQ        Frequency of operation, in Hz
%           FORMULATION Choice between 'JVIE_I' (1) or 'JVIE_II' (2)
%           LOCAL       if 0 uses the whole domain, otherwise only non-air
%           METHOD      Choice of an Iterative Solver
%                       0 if Self coded 'GMRES_DR' (with deflated restart)
%                       1 if Self coded 'GPU Dedicated JVIE Solver based on GMRES Deflated Restarted'
%                       2 if Self coded 'GMRES_R' (with restart)
%                       3 if Self coded 'GPU Dedicated JVIE Solver based on GMRES Restarted'
%                       4 if Self coded 'GMRES'
%                       5 if Matlab built-in 'M_GMRES'
%                       6 if Matlab built-in 'M_GMRES_R' (with restart)
%                       7 if Matlab built-in 'M_QMR'
%                       8 if Matlab built-in 'M_TFQMR'
%                       9 if Matlab built-in 'M_BICG'
%                       10 if Matlab built-in 'M_BICGSTAB'
%                       11 if Self coded 'GPU Dedicated BICGSTAB'
%           TOL         Relative tolerance for the method
%           IT          Number maximum iterations
%           INNER_IT    Number of internal iterations in restarted methods
%           OUTER_IT     Number of external iterations in restarted methods
%           RITZ        Number of vectors kept for deflated restart
%           GPU         if 0, no GPU, otherwise selects the GPU device
%           PRECOND     if 1, left preconditioning
%                       if 2, right preconditioning
%                       0 or otherwise, no preconditioner
%           LOGFILE     file to print the data
%                       
%
%   OUTPUT: JSOL        solution vector
%           FLAG        0 if converged, 1 if not
%           RELRES      final relative residue: norm(b - A*x)/norm(b) 
%           ITER        vector with [current internal iterations, current external iterations]
%           RESVEC      vector containing norm of relative residual at each iteration of GMRES
%           TIME        elapsed time by the method, in seconds
%
%
%
% _________________________________________________________________________
%
%   Computational Prototyping Group, RLE at MIT
% _________________________________________________________________________
% _________________________________________________________________________
%
%


% Initialize variables.
if(nargin < 5 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 6 || isempty(FORMULATION))
   FORMULATION = 2;
end
if(nargin < 7 || isempty(LOCAL))
   LOCAL = 1;
end
if(nargin < 8 || isempty(METHOD))
   METHOD = 0;
end
if(nargin < 9 || isempty(TOL))
   TOL = 1e-3;
else
    if(TOL < 1e-15) || (TOL >= 1)
        TOL = 1e-3;
    end
end
% % % These 4 are initialized below accordingly to method
% % if(nargin < 10 || isempty(IT))
% %    IT = 300;
% % end
% % if(nargin < 11 || isempty(INNER_IT))
% %    INNER_IT = 30;
% % end
% % if(nargin < 12 || isempty(OUTER_IT))
% %    OUTER_IT = 20;
% % end
% % if(nargin < 13 || isempty(RITZ))
% %    RITZ = 5;
% % end
if(nargin < 14 || isempty(GPU))
   GPU = 0;
end
if(nargin < 15 || isempty(PRECOND))
   PRECOND = 0;
end
if(nargin < 15 || isempty(LOGFILE))
   fid = 1;
else
    fid = fopen(LOGFILE, 'a');
    if (fid == -1)
        fid = 1; % standard output
    end
end

tini = tic;
fprintf(fid,'\n ------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------


% EM related
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*FREQ;
tau = 1j*omega*eo*(EPS-1);

% voxel volume
dx = R(2,1,1,1) - R(1,1,1,1);
dV = dx^3;

% Domain related
[L,M,N,~] = size(R);
nD = L*M*N; % number of variables in the system


if LOCAL
    idxS = find(abs(EPS(:))-1); %these are the indexes of the non-air positions
    idxS = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid
else
    idxS = 1:3*L*M*N;
end


% -------------------------------------------------------------------------
% Select the formulation
% -------------------------------------------------------------------------

if (FORMULATION == 1) % JVIE_I formulation

    % set the multiplier to convert fields to suitable input
    mult = [tau(:); tau(:); tau(:)]; % in 3D vector
    % get the input
    Vexc = mult.*EINC(:);
    
    if (GPU && gpuDeviceCount && (METHOD~=3) && (METHOD~=1) && (METHOD~=11)) % USE GPU ADAPTED FUNCTIONS
       
        dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
        reset(dev);
        
        % send data to GPU mem
        FG_gpu = gpuArray(FG);
        EPS_gpu = gpuArray(EPS);
                
        % define the main function
        fA   = @(J,transp_flag)mv_op_fft_Nop_1_gpu(J, FG_gpu, EPS_gpu, dV, transp_flag, idxS);
        % define the preconditioner
        fP   = @(J,transp_flag)mv_op_fft_Nop_1_gpu(J, FG_gpu, 1./EPS_gpu, dV, transp_flag, idxS);
    
    else
        
        % define the main function
        fA   = @(J,transp_flag)mv_op_fft_Nop_1(J, FG, EPS, dV, transp_flag, idxS);
        % define the preconditioner
        fP   = @(J,transp_flag)mv_op_fft_Nop_1(J, FG, 1./EPS, dV, transp_flag, idxS);
    
    end


else % default is the second formulation
    
    % set the multiplier to convert fields to suitable input
    mult = [tau(:)./EPS(:); tau(:)./EPS(:); tau(:)./EPS(:)]; % in 3D vector
    % get the input
    Vexc = mult.*EINC(:);
    
    if (GPU && gpuDeviceCount && (METHOD~=3) && (METHOD~=1) && (METHOD~=11)) % USE GPU ADAPTED FUNCTIONS
        
        dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
        reset(dev);
        
        % send data to GPU mem
        FG_gpu = gpuArray(FG);
        EPS_gpu = gpuArray(EPS);
                
        % define the main function
        fA   = @(J,transp_flag)mv_op_fft_Nop_2_gpu(J, FG_gpu, EPS_gpu, dV, transp_flag, idxS);
        % define the preconditioner
        fP   = @(J,transp_flag)mv_op_fft_Nop_2_gpu(J, FG_gpu, 1./EPS_gpu, dV, transp_flag, idxS);
        
    else
        
        % define the main function
        fA   = @(J,transp_flag)mv_op_fft_Nop_2(J, FG, EPS, dV, transp_flag, idxS);
        % define the preconditioner
        fP   = @(J,transp_flag)mv_op_fft_Nop_2(J, FG, 1./EPS, dV, transp_flag, idxS);
        
    end
    
end

Vexc = Vexc(idxS); % apply the selection


% -------------------------------------------------------------------------
% Select the iterative method and solve
% -------------------------------------------------------------------------
   

% -------------------------------------------------------------------------
% MATLAB BICGSTAB

if (METHOD==10) || (METHOD ==11) % 'M_BICGSTAB' or  'GPU JVIE BICGSTAB'

    if isempty(IT) || (IT < 1)
        initer = 300;
    else
        initer = IT;
    end
    
    if (PRECOND == 1)
        
        % MATLAB BICGSTAB: preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = bicgstab(@(J)fA(J,'notransp'), Vexc, TOL, initer, @(J)fP(J,'notransp'));
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN BICGSTAB FUNCTION, WITH PRECONDITIONER \n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));

       
    else
                
        if (METHOD == 11) && GPU && gpuDeviceCount
            
            % GPU Accelerated Dedicated JVIE solver Based on BICGSTAB
            
            tausolver = (EPS - 1.0) ./ EPS ; % value of tau
            dx3 = (R(2,1,1,1) - R(1,1,1,1))^3; % volume of the voxel
            
            dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
            reset(dev);
            
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = JJVIEsolverBicgstabGPU(FG, tausolver, L, M, N, dx3, idxS, Vexc, TOL, initer);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GPU ACCELERATED BICGSTAB FUNCTION, NO PRECONDITIONER \n');
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, relative residue %g \n' ,method_time,flag_,length(resvec_),resvec_(end));
            
        else            
            
            % MATLAB BICGSTAB: no preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = bicgstab(@(J)fA(J,'notransp'), Vexc, TOL, initer);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n MATLAB BUILT IN BICGSTAB FUNCTION, NO PRECONDITIONER\n');
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        end

        
    end
        
end



% -------------------------------------------------------------------------
% MATLAB BICG

if (METHOD==9) % 'M_BICG'

    if isempty(IT) || (IT < 1)
        initer = 300;
    else
        initer = IT;
    end
    
    if (PRECOND == 1)
        
        % MATLAB BICG: preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = bicg(fA, Vexc, TOL, initer, fP);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN BICG FUNCTION, WITH PRECONDITIONER \n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
              
    else
        
        % MATLAB BICG: no preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = bicg(fA, Vexc, TOL, initer);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN BICG FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
        
    end
        
end



% -------------------------------------------------------------------------
% MATLAB TFQMR

if (METHOD==8) % 'M_TFQMR'

    if isempty(IT) || (IT < 1)
        initer = 300;
    else
        initer = IT;
    end
    
    if (PRECOND == 1)
        
        % MATLAB TFQMR: preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = tfqmr(@(J)fA(J,'notransp'), Vexc, TOL, initer, @(J)fP(J,'notransp'));
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN TFQMR FUNCTION, WITH PRECONDITIONER \n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));

    else            
        
        % MATLAB TFQMR: no preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = tfqmr(@(J)fA(J,'notransp'), Vexc, TOL, initer);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN TFQMR FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
  
    end
        
end



% -------------------------------------------------------------------------
% MATLAB QMR

if (METHOD==7)

    if isempty(IT) || (IT < 1)
        initer = 300;
    else
        initer = IT;
    end
    
    if (PRECOND == 1)
        
        % MATLAB QMR: preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = qmr(fA, Vexc, TOL, initer, fP);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN QMR FUNCTION, WITH PRECONDITIONER \n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));

    else
                        
        % MATLAB QMR: no preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = qmr(fA, Vexc, TOL, initer);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN QMR FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
        
    end
        
end



% -------------------------------------------------------------------------
% MATLAB GMRES

if (METHOD==5) || (METHOD==6) 

    if (METHOD==5) % NO RESTARTED VERSION
        if isempty(IT) || (IT < 1)
            initer = 300;
        else
            initer = IT;
        end        
        outiter = 1;
    else % RESTARTED METHOD
        if isempty(INNER_IT) || (INNER_IT < 1)
            initer = 30;
        else
            initer = INNER_IT;
        end        
        if isempty(OUTER_IT) || (OUTER_IT < 1)
            outiter = 40;
        else
            outiter = OUTER_IT;
        end
    end
    
    if (PRECOND == 1)
        
        % GMRES: left preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = gmres(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, @(J)fP(J,'notransp'));
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN GMRES (%d) FUNCTION, LEFT PRECONDITIONER \n', initer );
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
   
    else
        
        % GMRES: no preconditioning
        ta = tic;
        [vsol_,flag_,relres_,iter_,resvec_] = gmres(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter);
        method_time = toc(ta);
        r_ = Vexc - fA(vsol_ , 'notransp');
        fprintf(fid,'\n MATLAB BUILT IN GMRES (%d) FUNCTION, NO PRECONDITIONER \n', initer );
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
                
    end
        
end



% -------------------------------------------------------------------------
% SELF IMPLEMENTATION OF GMRES

if  (METHOD == 2) || (METHOD == 3) || (METHOD == 4)

    if (METHOD==4) % NO RESTARTED VERSION
        if isempty(IT) || (IT < 1)
            initer = 300;
        else
            initer = IT;
        end        
        outiter = 1;
    else % RESTARTED METHOD
        if isempty(INNER_IT) || (INNER_IT < 1)
            initer = 30;
        else
            initer = INNER_IT;
        end        
        if isempty(OUTER_IT) || (OUTER_IT < 1)
            outiter = 40;
        else
            outiter = OUTER_IT;
        end
    end
    
    if (PRECOND < 1) || (PRECOND > 2) % NO PRECONDITIONING
            
        if (METHOD == 3) && GPU && gpuDeviceCount
            
            % GPU Accelerated Dedicated JVIE solver Based on GMRES 
            
            tausolver = (EPS - 1.0) ./ EPS ; % value of tau
            dx3 = (R(2,1,1,1) - R(1,1,1,1))^3; % volume of the voxel
            
            dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
            reset(dev);
            
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = JVIEsolverGPU(FG, tausolver, L, M, N, dx3, idxS, Vexc, initer, TOL, outiter);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GPU ACCELERATED GMRES (%d) FUNCTION, NO PRECONDITIONER \n', initer );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, relative residue %g \n' ,method_time,flag_,length(resvec_),resvec_(end));
            
        else
            
            % GMRES: no preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmres(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES (%d) FUNCTION, NO PRECONDITIONER \n', initer );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        end
        
    else
        
        if (PRECOND < 2)
            
            % GMRES: left preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmres(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, @(J)fP(J,'notransp'));
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES (%d) FUNCTION, LEFT PRECONDITIONER \n', initer );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        else
            
            % GMRES: right preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmres(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, [], [], @(J)fP(J,'notransp'));
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES (%d) FUNCTION, RIGHT PRECONDITIONER \n', initer );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        end
                
    end
        
end


% -------------------------------------------------------------------------
% GMRES DR

if (METHOD < 2) || (isempty(vsol_)) % DEFAULT

    if isempty(INNER_IT) || (INNER_IT < 1)
        initer = 30;
    else
        initer = INNER_IT;
    end
    
    if isempty(OUTER_IT) || (OUTER_IT < 1)
        outiter = 40;
    else
        outiter = OUTER_IT;
    end
    
    if isempty(RITZ) || (RITZ < 1)
        RITZ = 5;
    else
        RITZ = min(RITZ,initer);
    end
    
    
    if (PRECOND < 1) || (PRECOND > 2) % NO PRECONDITIONING
            
        if (METHOD == 1) && GPU && gpuDeviceCount
            
            % GPU accelerated dedicated JVIE solver Based on GMRES-DR
            
            tausolver = (EPS - 1.0) ./ EPS ; % value of tau
            dx3 = (R(2,1,1,1) - R(1,1,1,1))^3; % volume of the voxel
            
            dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
            reset(dev);

            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = JVIEsolverDRGPU(FG, tausolver, L, M, N, dx3, idxS, Vexc, initer, TOL, outiter, RITZ);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GPU ACCELERATED GMRES-DR (%d,%d) FUNCTION, NO PRECONDITIONER \n', initer, RITZ );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        else
                        
            % GMRES DR: no preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmresDR(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, RITZ);
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES DR (%d,%d) FUNCTION, NO PRECONDITIONER \n', initer, RITZ );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        end
        
    else
        
        if (PRECOND < 2)
            
            % GMRES DR: left preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmresDR(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, RITZ, @(J)fP(J,'notransp'));
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES DR (%d,%d) FUNCTION, LEFT PRECONDITIONER \n', initer, RITZ );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        else
            
            % GMRES DR: right preconditioning
            ta = tic;
            [vsol_,flag_,relres_,iter_,resvec_] = pgmresDR(@(J)fA(J,'notransp'), Vexc, initer, TOL, outiter, RITZ, [], [], @(J)fP(J,'notransp'));
            method_time = toc(ta);
            r_ = Vexc - fA(vsol_ , 'notransp');
            fprintf(fid,'\n GMRES DR (%d,%d) FUNCTION, RIGHT PRECONDITIONER \n', initer, RITZ );
            fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %f, relative residue %g \n' ,method_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(Vexc));
            
        end
        
    end

end



% -------------------------------------------------------------------------
% Assign veriables for return and clean whatever created... in GPU
% -------------------------------------------------------------------------


if (GPU && gpuDeviceCount) % USE GPU ADAPTED FUNCTIONS
    
    dev = gpuDevice(GPU); % assign the GPU selected (1, 2, ... up to gpuDeviceCount)
    reset(dev);
    
end


JSOL = zeros(3*L*M*N,1);
JSOL(idxS) = vsol_ ; % return to global variables
RELRES = relres_ ;
ITER = iter_ ;
RESVEC = resvec_;
TIME = method_time;

tend = toc(tini);

if (norm(r_)/norm(Vexc)) > TOL
    FLAG = 1; % no real convergence
    fprintf(fid,'\n SOLVE FINISHED (warning: convergence criteria was not achieved). Elapsed Time  = %g [sec] \n' ,tend);
else
    FLAG = 0;
    fprintf(fid,'\n SOLVE FINISHED. Elapsed Time  = %g [sec] \n' ,tend);
end

if (fid ~= 1)
    fclose(fid);
end





