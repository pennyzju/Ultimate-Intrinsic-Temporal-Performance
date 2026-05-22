function [M] = MRGF_solution(M, e_r, r, f, outfile, tol)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates the MRGF for a given domain and eps
%
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%

if(nargin < 6)
    tol = 1e-4;
end

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

fid = 1;
% the output file will be in the ..\DATA_MRGF folder
outpath = strcat('.\DATA_MRGF\', outfile);

% -------------------------------------------------------------------------
%                  EM constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * f;
% lambda  = co/f;
% ko = 2*pi/lambda;
% omega_mu = omega*mu;
% eta =  3.767303134617706e+002; % Free-space impedance


% -------------------------------------------------------------------------
%            DATA
% -------------------------------------------------------------------------

% idxS = find(abs(e_r(:)-1));
dx = r(2,1,1,1) - r(1,1,1,1);
dV = dx^3;

% -------------------------------------------------------------------------
%                  Prepare the cases to Solve
% -------------------------------------------------------------------------

% % define subdomain of interest
% xd = r(:,:,:,1);
% xmin = min(xd(idxS))-dx; xmax = max(xd(idxS))+dx;
% yd = r(:,:,:,2);
% ymin = min(yd(idxS))-dx; ymax = max(yd(idxS))+dx;
% zd = r(:,:,:,3);
% zmin = min(zd(idxS))-dx; zmax = max(zd(idxS))+dx;
% 
% % cut grid to subdomain dimensions
% idxX = find( (r(:,1,1,1) >= xmin) & (r(:,1,1,1) <= xmax));
% idxY = find( (r(1,:,1,2) >= ymin) & (r(1,:,1,2) <= ymax));
% idxZ = find( (r(1,1,:,3) >= zmin) & (r(1,1,:,3) <= zmax));


% -------------------------------------------------------------------------
%                  Create new grid
% -------------------------------------------------------------------------

% % truncate e_r to new reduced grid
% e_rnew = e_r(idxX,idxY,idxZ);
% 
e_rnew = e_r;

% get dimensions
[Lnew,Mnew,Nnew] = size(e_rnew);

% -------------------------------------------------------------------------
%                  Determine good circulant size and compute if needed
% -------------------------------------------------------------------------

% Multiples of 2
LfG = 2*(ceil(Lnew/2));
MfG = 2*(ceil(Mnew/2));
NfG = 2*(ceil(Nnew/2));

% get the required circulant for the EF operator
[fG] = retrieve_fG(LfG,MfG,NfG,dx,f);

% -------------------------------------------------------------------------
%                 Define the solution parameters
% -------------------------------------------------------------------------

% prepare domain
idxE = find( abs(e_rnew(:))-1 );
idxE3 = [idxE; Lnew*Mnew*Nnew+idxE; 2*Lnew*Mnew*Nnew+idxE]; % the vector of non-air positions in 3D grid


infocir = whos('fG');

if (gpuDeviceCount)
    dev = gpuDevice;
    
    if (dev.TotalMemory > 20*infocir.bytes/6)
        reset(dev);

        % send data to GPU mem
        FG_gpu = gpuArray(fG);
        EPS_gpu = gpuArray(e_rnew);
                
        % define the main function
        fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying GPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fA   = @(J,sTrans)mv_op_fft_Nop_2_gpu(J, FG_gpu, EPS_gpu, dV, sTrans, idxE3);
        
    else        
        fprintf(1,'\n Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
        fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_rnew, dV, sTrans, idxE3);
        
    end        
    
end
    
if ~(exist('fA'))
    
    fprintf(1,'\n No GPU available. Estimated MVP Peak Memory %.3f MB.\n Applying CPU based Solver \n' , 12*infocir.bytes/(6*1024*1024));
    fA   = @(J,sTrans)mv_op_fft_Nop_2(J, fG, e_rnew, dV, sTrans, idxE3);
    
end



% -------------------------------------------------------------------------
% Obtain the Scattered field due to Qe
% -------------------------------------------------------------------------

% tol = 1e-4;
maxit = 1000; % number of external iterations

% set the multiplier to convert fields to suitable input
tau = 1j*omega*eo*(e_rnew(:)-1)./(e_rnew(:));
mult = [tau; tau; tau];
mult = dV*mult(idxE3); % chose Scatterer positions

% allocate space
ncol = size(M,2);

fprintf(fid, '\n  Solving the system for %d inputs \n', ncol);
t1 = tic;

for ii=1:ncol

    % get the incident field excitation and transform into currents
    JExc = mult.*M(:,ii);
         
    % solve to obtain the total currents
    [jsol,~,~,~,~] = bicgstab(@(J)fA(J,'notransp'), JExc, tol, maxit);
    
    M(:,ii) = jsol; % store the solution currents
                  
    if ((ii/500) == floor(ii/500))
        fprintf(fid, '\n %d Solves. Elapsed time %g', ii, toc(t1));
    end
        
    
end

clear Jout; clear jsol; clear JExc;

fprintf(fid, '\n\n');

t2 = toc(t1);
fprintf(fid, '\n\n Equivalent current generated for %d vectors. Elapsed time %g \n', ncol, t2);

Xout_file = sprintf('%s_MRGF_M.mat', outpath);
save(Xout_file, 'M', '-v7.3');

