function [Jout,Sout,Eout,Bout,fG,fK] = solver_all(Ein, Hin, f, r, Eps_r, Sig_e, Rho, out, tol, bicgtrue, nogpu)
% _________________________________________________________________________
% _________________________________________________________________________
%
%  JVIE method for solving the system
%
% _________________________________________________________________________
% _________________________________________________________________________


if(nargin < 7)
   Rho = [];
end

if(nargin < 8)
   out = 1;
else
    if out < 1
        out = 1;
    end        
end

if(nargin < 9)
   tol = 1e-4;
end

if(nargin < 10)
   bicgtrue = 1;
end

if(nargin < 11)
   nogpu = 0;
end

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * f;
lambda  = co/f;
ko = 2*pi/lambda;

dx = r(2,1,1,1) - r(1,1,1,1);
dV = dx^3;

Tesla = f/(42.6e6);

e_r = Eps_r - (1j/(eo * omega)) * Sig_e;


% -------------------------------------------------------------------------
%                  Prepare the cases to Solve
% -------------------------------------------------------------------------

[L,M,N] = size(e_r); % get dimensions of the problem
nD = L*M*N; % number of voxels in the complete domain

[~, maxeps] = max(real(e_r(:)));
idxS = find(abs(e_r(:)-1));
idxSS = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid

% Convert the point fields into equivalent volumetric values
Ein = Ein*dV;
Hin = Hin*dV;

% translate the input in the dimensions of the scatterer
% set the multiplier to convert fields to suitable input
tau = 1j*omega*eo*(e_r-1)./(e_r);
% get the input
Vexc = zeros(length(idxS),3);
for ii = 1:3
    Vcomp = tau.*Ein(:,:,:,ii);
    Vexc(:,ii) = Vcomp(idxS);
end
Vexc = Vexc(:);
clear Vcomp; clear tau;


fid = 1;
fprintf(fid, '\n\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n     # DOFS:                %d', 3*nD);
fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
fprintf(fid, '\n     Operating Frequency:   %.2f MHz (%.1f T)',f/1e6, Tesla);
fprintf(fid, '\n     Maximum Epsilon:       %.2f %.2f j', real(e_r(maxeps)), imag(e_r(maxeps)));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');


% -------------------------------------------------------------------------
%            Subdomain
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

% get dimensions
[L,M,N] = size(e_rnew);

% -------------------------------------------------------------------------
%                  Determine good circulant size and compute if needed
% -------------------------------------------------------------------------

% Multiples of 4
LfG = 4*(ceil(L/4));
MfG = 4*(ceil(M/4));
NfG = 4*(ceil(N/4));

% get the required circulant
[fG] = retrieve_fG(LfG,MfG,NfG,dx,f);

% -------------------------------------------------------------------------
%                 Solve the system for the croped domain
% -------------------------------------------------------------------------


% prepare domain
idxE = find(abs(e_rnew(:))-1);
idxEE = [idxE; L*M*N+idxE; 2*L*M*N+idxE]; % the vector of non-air positions in 3D grid


% BICGSTAB: no preconditioning
maxit = 5000;
tic_solve = tic;

% function that computes the scattered field due to an excitation field, including
% non-air media
infocir = whos('fG');

if (gpuDeviceCount) && (nogpu == 0)
    dev = gpuDevice;
    
    if (dev.TotalMemory > 20*infocir.bytes/6)
        reset(dev);

        % send data to GPU mem
        FG_gpu = gpuArray(fG);
        EPS_gpu = gpuArray(e_rnew);
                
        % define the main function
        fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fA   = @(J,transp_flag)mv_op_fft_Nop_2_gpu(J, FG_gpu, EPS_gpu, dV, transp_flag, idxEE);
        
    else        
        
        if (dev.TotalMemory > 7*infocir.bytes/6)
            reset(dev);
            
            fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying semi-GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
            fA   = @(J,sTrans)mv_op_fft_Nop_2_semigpu(J, fG, e_rnew, dV, sTrans, idxEE);
            
        else
                        
            fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
            fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_rnew, dV, sTrans, idxEE);
            
        end
        
    end        
    
end
    
if ~(exist('fA'))
    
    fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n No GPU, Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
    fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_rnew, dV, sTrans, idxEE);
    
end

if bicgtrue
    
    [vSol,flag,~,~,resvec] = bicgstab(@(J)fA(J,'notransp'), Vexc, tol, maxit);
    % [vSol,flag,~,~,resvec] = Mbicgstab(@(J)fA(J,'notransp'), Vexc, tol, maxit,[],idxEE,L,M,N,round(L/2),round(M/2),3*round(N/4));
    
    Time= toc(tic_solve);
    residue =  Vexc - fA(vSol,'notransp');
    fprintf(1,'\n Current Distribution Computed:\n BICGSTAB Elapsed time  = %g [sec], flag %d, %d iterations, residue %f, relative %g \n' , Time,flag,length(resvec),norm(residue),norm(residue)/norm(Vexc));
        
else
    
    restart = 50;
    [vSol,flag,~,~,resvec] = pgmres(@(J)fA(J,'notransp'), Vexc, restart, tol, ceil(maxit/restart));
    
    Time= toc(tic_solve);
    residue =  Vexc - fA(vSol,'notransp');
    fprintf(1,'\n Current Distribution Computed:\n GMRES Elapsed time  = %g [sec], flag %d, %d iterations, residue %f, relative %g \n' , Time,flag,length(resvec),norm(residue),norm(residue)/norm(Vexc));
    
end


clear e_rnew; clear idxE; clear idxEE;

% -------------------------------------------------------------------------
%                 Translate to the complete original domain
% -------------------------------------------------------------------------

[L,M,N] = size(e_r);

% translate solution to the complete domain
Jout = zeros(L,M,N,3);
Jout(idxSS) = vSol;
Jout = reshape(Jout,L,M,N,3);

% -------------------------------------------------------------------------
%                   Compute SAR
% -------------------------------------------------------------------------

if (out > 1)
    
    tic_rad = tic;
    
    tau = 1j*omega*eo*(e_r-1);
    Sout = zeros(L,M,N);
    
    for ii = 1:3
        Jcomp = Jout(:,:,:,ii);
        Sout(idxS) = Sout(idxS) + (Jcomp(idxS)./tau(idxS)).*conj(Jcomp(idxS)./tau(idxS));
    end
    
    Sout(idxS) = Sout(idxS).*Sig_e(idxS)./(2*Rho(idxS));
    
    Time= toc(tic_rad);
    fprintf(1,'\n Voxelized SAR computed:\n  Elapsed time  = %g [sec]\n' , Time);
    
else
    Sout = [];
end
    

% -------------------------------------------------------------------------
%                   Compute E field
% -------------------------------------------------------------------------

if (out > 2)
    
    % get the required circulant
    [fG] = retrieve_fG(L,M,N,dx,f);
    
    tic_rad = tic;
    
    % generate the radiated field from the currents
    [Eout] = E_field_Nop_comp(Jout, fG, dV, omega, eo, Ein);
    
    Time= toc(tic_rad);
    fprintf(1,'\n Electric Field computed:\n  Elapsed time  = %g [sec]\n' , Time);
    
    % add the incident field to the radiated by the currents
    Eout = reshape(Eout,L,M,N,3);
    
else
    Eout = [];
end


% -------------------------------------------------------------------------
%                   Compute H field
% -------------------------------------------------------------------------

if (out > 3)
        
    % get the required circulant
    [fK] = retrieve_fK(L,M,N,dx,f);

    tic_rad = tic;
    
    % generate the B field from the currents
    [Hout] = H_field_Kop_comp(Jout, fK, dV, Hin);
    Bout = mu*Hout;
    
    Time= toc(tic_rad);
    fprintf(1,'\n Magnetic Field computed:\n  Elapsed time  = %g [sec]\n' , Time);
    
    % add the incident field to the radiated by the currents
    Bout = reshape(Bout,L,M,N,3);
   
    
else
    Bout = [];
    fK = [];
end





% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');