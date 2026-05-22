function [Eout,Jout,fG] = jvie_solverbis(Vin, e_r, r, res, f, fG, tol, x0)
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
if(nargin < 8)
   x0 = [];
end

% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------


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
% fprintf(fid, '\n\n ----------------------------------------------------------');
% fprintf(fid, '\n ----------------------------------------------------------');
% fprintf(fid, '\n');
% fprintf(fid, '\n');
% fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
% fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
% fprintf(fid, '\n     # DOFS:                %d', 3*nD);
% fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
% fprintf(fid, '\n     Operating Frequency:   %.2f MHz (%.1f T)',f/1e6,f/(42.6e6));
% fprintf(fid, '\n     Maximum Epsilon:       %.2f %.2f j', real(e_r(maxeps)), imag(e_r(maxeps)));
% fprintf(fid, '\n');
% fprintf(fid, '\n ----------------------------------------------------------');
% fprintf(fid, '\n ----------------------------------------------------------\n\n ');


% -------------------------------------------------------------------------
%                  Compute the circulant for the croped domain
% -------------------------------------------------------------------------

if isempty(fG) % if there is no Gmn, create it

    t1 = clock;
    
    % Generate Green's Tensor G_mn
    [Gmn] = Generate_Gmn(r,ko);
    
    t2 = clock;
    fprintf(fid, '\n Green Tensor computed. elapsed time %f', etime(t2,t1));
    
    t1 = clock;
    
    % Generate Circulant matrix
    [fG] = Generate_Circulant(Gmn);
    
    fprintf(1,'\n Circulant computed. Elapsed time %g [sec]\n', etime(clock,t1));
    
else % get the dimensions we need
    
%     fprintf(1,'\n Circulant Loaded\n');
    
end


% save('temp_fG.mat', 'fG', '-V7.3');

% infocir = whos('fG');
% fprintf(fid, '\n ----------------------------------------------------------');
% fprintf(fid, '\n     FFT Space Dimensions:   %dx%dx%d',size(fG,1),size(fG,2), size(fG,3));
% fprintf(fid, '\n     Circulant Memory:       %.6f MB', infocir.bytes/(1024*1024));
% fprintf(fid, '\n     Estimated Peak Memory:  %.6f MB', 12*infocir.bytes/(6*1024*1024));
% fprintf(fid, '\n');
% fprintf(fid, '\n ----------------------------------------------------------\n\n ');


% -------------------------------------------------------------------------
%                 Solve the system for the croped domain
% -------------------------------------------------------------------------


% prepare domain
idxSS = [idxS; L*M*N+idxS; 2*L*M*N+idxS]; % the vector of non-air positions in 3D grid

infocir = whos('fG');

if (gpuDeviceCount)
    dev = gpuDevice;
    
    if (dev.TotalMemory > 20*infocir.bytes/6)
        reset(dev);

        % send data to GPU mem
        FG_gpu = gpuArray(fG);
        EPS_gpu = gpuArray(e_r);
                
        % define the main function
%         fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fA   = @(J,sTrans)mv_op_fft_Nop_2_gpu(J, FG_gpu, EPS_gpu, dV, sTrans, idxSS);
        
    else        
%         fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_r, dV, sTrans, idxSS);
        
    end        
    
end
    
if ~(exist('fA'))
    
%     fprintf(1,'\n No GPU available. Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
    fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_r, dV, sTrans, idxSS);
    
end

% set the multiplier to convert fields to suitable input
tau = 1j*omega*eo*(e_r-1)./(e_r);
% get the input
Vexc = zeros(L,M,N,3);
for ii = 1:3
    Vexc(:,:,:,ii) = tau.*Vin(:,:,:,ii);
end

% fprintf(1,'\n Solving with the JVIE II formulation\n');

% GMRES DR: no preconditioning
maxit = 1000;
tic_bicgstab = tic;

[vSol,flag,~,~,resvec] = bicgstab(@(J)fA(J,'notransp'), Vexc(idxSS), tol, maxit);
%[vSol,flag,~,~,resvec] = Mbicgstab(@(J)fA(J,'notransp'), Vexc(idxSS), tol, maxit,x0,idxSS,L,M,N,round(L/2),round(M/2),3*round(N/4));

Time= toc(tic_bicgstab);
residue =  Vexc(idxSS) - fA(vSol,'notransp');
fprintf(1,'\n Current Distribution Computed:\n  Elapsed time  = %g [sec], flag %d, %d iterations, residue %f, relative %g \n' , Time,flag,length(resvec),norm(residue),norm(residue)/norm(Vexc(idxSS)));

  

% -------------------------------------------------------------------------
%                 Translate to the complete original domain
% -------------------------------------------------------------------------


% translate solution to the complete domain
Jout = zeros(L,M,N,3);
Jout(idxSS) = vSol;

tic_rad = tic;

% generate the radiated field from the currents
[Eout] = field_comp(Jout, fG, dV, omega, eo, Vin);

Time= toc(tic_rad);
% fprintf(1,' Field computed:  Elapsed time  = %g [sec]\n' , Time);

% add the incident field to the radiated by the currents
Eout = reshape(Eout,L,M,N,3);
Jout = reshape(Jout,L,M,N,3);

% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

% fprintf(fid, '\n');
% fprintf(fid, '\n ----------------------------------------------------------');
% fprintf(fid, '\n ----------------------------------------------------------\n\n ');