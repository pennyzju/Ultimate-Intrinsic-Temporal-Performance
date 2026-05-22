function [Eout,Jout,fG] = jvie_solver_v2(Vin, e_r, r, res, f, fG, tol)
% _________________________________________________________________________
% _________________________________________________________________________
%
%  JVIE method for solving the system
%
% _________________________________________________________________________
% _________________________________________________________________________


if(nargin < 6)
   fG = []; % empty fG... to compute the circulant
end
if(nargin < 7)
   tol = 1e-3;
end

% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------

global ko
global dx 
global R_faces   %#ok<NUSED>

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * f;
lambda  = co/f;
ko = 2*pi/lambda;

dx = res;
dV = dx^3;


% -------------------------------------------------------------------------
%                  Prepare the cases to Solve
% -------------------------------------------------------------------------

[L,M,N] = size(e_r); % get dimensions of the problem
nD = L*M*N; % number of voxels in the complete domain

[~, maxeps] = max(real(e_r(:)));
idxS = find(abs(e_r(:))-1);

fid = 1;
fprintf(fid, '\n\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n     # DOFS:                %d', 3*nD);
fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
fprintf(fid, '\n     Operating Frequency:   %.2f MHz (%.1f T)',f/1e6,f/(42.6e6));
fprintf(fid, '\n     Maximum Epsilon:       %.2f %.2f j', real(e_r(maxeps)), imag(e_r(maxeps)));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');


% -------------------------------------------------------------------------
%                  Compute the circulant for the croped domain
% -------------------------------------------------------------------------

if isempty(fG) % if there is no Gmn, create it

    t1 = clock;
    
    % Generate Green's Tensor G_mn
    [Gmn] = Generate_Gmn(r);
    
    t2 = clock;
    fprintf(fid, '\n Green Tensor computed. elapsed time %f', etime(t2,t1));
    
    t1 = clock;
    
    % Generate Circulant matrix
    [fG] = Generate_Circulant(Gmn);
    
    fprintf(1,'\n Circulant computed. Elapsed time %g [sec]\n', etime(clock,t1));
    
else % get the dimensions we need
    
    fprintf(1,'\n Circulant Loaded\n');
    
end


% save('temp_fG.mat', 'fG', '-V7.3');

infocir = whos('fG');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n     FFT Space Dimensions:   %dx%dx%d',size(fG,1),size(fG,2), size(fG,3));
fprintf(fid, '\n     Circulant Memory:       %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n     Estimated Peak Memory:  %.6f MB', 12*infocir.bytes/(6*1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');


% -------------------------------------------------------------------------
%                 Solve the system for the croped domain
% -------------------------------------------------------------------------


% prepare domain
% idxS = [idxS; L*M*N+idxS; 2*L*M*N+idxS]; % the vector of non-air positions in 3D grid
idxS = 1:3*L*M*N;
idxS = idxS.';


% set the multiplier to convert fields to suitable input
tau = 1j*omega*eo*(e_r-1)./(e_r);
% get the input
Vexc = zeros(L,M,N,3);
for ii = 1:3
    Vexc(:,:,:,ii) = tau.*Vin(:,:,:,ii);
end


notsolved = 1;
% GMRES DR: no preconditioning
ritz = 5; maxit = 50; maxout = 50;
tic_drgmres = tic;

if (0)
    dev = gpuDevice;
    
    if (dev.TotalMemory > 12*infocir.bytes/6)
        reset(dev);
        
        % tausol = (e_rnew-1)./(e_rnew);        
        % [vSol,flag,~,~,resvec] = JVIEsolverGPU(fG, tausol, L, M, N, dV, idxE, Vexc, maxit, tol, maxout);
                
        % send data to GPU mem
        FG_gpu = gpuArray(fG);
        EPS_gpu = gpuArray(e_r);
        
        fprintf(1,'\n Solving on GPU with the JVIE II formulation\n');
        tic_drgmres = tic;
        % define the main function
        fA   = @(J,transp_flag)mv_op_fft_2_gpu(J, FG_gpu, EPS_gpu, dV, transp_flag, idxS);
        [vSol,flag,~,~,resvec] = bicgstab(@(J)fA(J,'notransp'), Vexc(:), tol, maxit*maxout);
        
        Time= toc(tic_drgmres);
        residue =  Vexc - fA(vSol,'notransp');
        fprintf(1,' GPU Solver:  Elapsed time  = %g [sec], flag %d, %d iterations, residue %f, relative %g \n' , Time,flag,length(resvec),norm(residue),norm(residue)/norm(Vexc));
        
        notsolved = 0;
        
    end
    
end
    
if (notsolved)
    % function that computes the scattered field due to an excitation field, including
    % non-air media
    
    fprintf(1,'\n Solving with the JVIE II formulation\n');
    fA  = @(J,sTrans)mv_op_fft_2(J, fG, e_r, dV, sTrans, idxS);
    [vSol,flag,~,~,resvec] = pgmresDR(@(J)fA(J,'notransp'), Vexc(:), maxit, tol, maxout, ritz);
    
    Time= toc(tic_drgmres);
    residue =  Vexc - fA(vSol,'notransp');
    fprintf(1,' Solver:  Elapsed time  = %g [sec], flag %d, %d iterations, residue %f, relative %g \n' , Time,flag,length(resvec),norm(residue),norm(residue)/norm(Vexc));
    
end

  

% -------------------------------------------------------------------------
%                 Translate to the complete original domain
% -------------------------------------------------------------------------


% translate solution to the complete domain
Jout = zeros(L,M,N,3);
Jout(idxS) = vSol;

tic_rad = tic;

% generate the radiated field from the currents
[Eout] = field_comp(Jout, fG, dV, omega, eo, Vin);

Time= toc(tic_rad);
fprintf(1,' Field computed:  Elapsed time  = %g [sec]\n' , Time);

% add the incident field to the radiated by the currents
Eout = reshape(Eout,L,M,N,3);
Jout = reshape(Jout,L,M,N,3);

% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');