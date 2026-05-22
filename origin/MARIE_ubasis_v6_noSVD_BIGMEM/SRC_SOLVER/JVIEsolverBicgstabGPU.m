% GPU adapted Dedicated Solver for JVIE formulation
% using BICGSTAB, with halfstep convergence checking
%
%   Version by J. Fernandez Villena
%              Computational Prototyping Group, RLE at MIT
%
% INPUT:  fG        Circulant matrix for function application
%         tau       inhomogeneous properties at given frequency
%         L         number of discretization elements in x direction
%         M         number of discretization elements in y direction
%         N         number of discretization elements in z direction
%         dx3       volume of discretized voxel
%         idx       indexes of the elements relevant for the solution
%         b         rhs
%         tol       relative tolerance for the residue
%         maxit     maximum number of outer iterations (total is restart*maxit)
%         x         initial guess
%
% OUTPUT: x         solution vector
%         flag      0 if converged, 1 if not
%         relres    final relative residue: norm(b - A*x)/norm(b) 
%         iter      vector with [current iteration, solution iteration]
%         resvec    vector containing norm of relative residual at each half iteration
%         
%
%
% SUBFUNCTIONS: 
%
%



function [x,flag,relres,iter,resvec] = JVIEsolverBicgstabGPU(fG, tau, L, M, N, dx3, idx, b, tol, maxit, x)

% Initialize variables.
if(nargin < 8 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 9 || isempty(tol))
   tol = 1e-3;
end
if(nargin < 10 || isempty(maxit))
   maxit = 500;
end
if(nargin < 11 || isempty(x))
    x = zeros(size(b));
end


% -------------------------------------------------------------------------
% settings

% send data to GPU
fG = gpuArray(fG);
tau = gpuArray(tau);
idx = gpuArray(idx);
b = gpuArray(b);
x = gpuArray(x);


% Calculate rhs norm
bnorm = norm(b);

% initialize the residue vector 
resvec = gpuArray.zeros(2*maxit+1,1);

% initialize function data
JIn0 = gpuArray.zeros(L, M, N, 3);
JOut = gpuArray.zeros(L, M, N, 3);

% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);


% -------------------------------------------------------------------------
% compute the initial residual
%   Apply the function A*x
%
JIn0(idx) = x;

% apply fft and mv-op
fJ = fftn(JIn0(:,:,:,1),[LfG, MfG, NfG]);
Jout1 = fG(:,:,:,1) .* fJ;
Jout2 = fG(:,:,:,2) .* fJ;
Jout3 = fG(:,:,:,3) .* fJ;

fJ = fftn(JIn0(:,:,:,2),[LfG, MfG, NfG]);
Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
Jout3 = Jout3 + fG(:,:,:,5) .* fJ;

fJ = fftn(JIn0(:,:,:,3),[LfG, MfG, NfG]);
Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
Jout3 = Jout3 + fG(:,:,:,6) .* fJ;

% apply ifft
Jout1 = ifftn(Jout1);
JOut(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
Jout2 = ifftn(Jout2);
JOut(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
Jout3 = ifftn(Jout3);
JOut(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);

% sum contribution and return local indexes
r = b - (dx3 * x - JOut(idx));

% -------------------------------------------------------------------------

% initialize vars for the BICGSTAB method

p = r; % initial p
ro = r; % initial residual
rho = ro'*r; % % compute norm of current residue
stag = 0; % stagnation flag
flag = 1; % convergence flag

bestii = 1; % best case residual norm
bestsol = x; % best case solution
resvec(1) = norm(r)/bnorm; % residual norm

for ii = 0:maxit

    % ---------------------------------------------------------------------
    % first half iteration
    % ---------------------------------------------------------------------
        
    % ---------------------------------------------------------------------
    %   Apply the function v = A*p on previous p 
    %
    
    JIn0(idx) = p;
    
    % apply fft and mv-op
    fJ = fftn(JIn0(:,:,:,1),[LfG, MfG, NfG]);
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn0(:,:,:,2),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn0(:,:,:,3),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);
    
    % sum contribution and return local indexes
    v = dx3 * p - JOut(idx);
    % v = A*p; % apply the operation    
    % ---------------------------------------------------------------------
    
        
    alpha = rho / (ro'*v); % compute alpha coeff
    x = x +  alpha*p; % update solution
    s = r - alpha*v; % generate s (shadow residual)
    currentii = 2*(ii+1);
    resvec(currentii) = norm(s)/bnorm; % store half iteration residual
    
    % check if better solution that current one
    if resvec(currentii) < resvec(bestii)
        
        bestii = currentii;
        
        % if better solution check half convergence
        if ( resvec(bestii) <= tol)
            flag = 0;
            break;
        end
        
        % if not, store as best solution so far
        bestsol = x;
        stag = 0;
        
    else
        
        % check for stagnation or bad conditioning of variables
        if abs(alpha) < 1e-12
            if (stag > 5)
                flag = 3; % stagnated
                break
            else
                stag = stag + 1;
            end
        else
            if isfinite(alpha)
                stag = 0;
            else
                flag = 4; % numerical errors
                break
            end
        end

        
    end
    
    
    % ---------------------------------------------------------------------
    % Second half iteration
    % ---------------------------------------------------------------------
        
    % ---------------------------------------------------------------------
    %   Apply the function t = A*s on shadow residual
    %
    
    JIn0(idx) = s;
    
    % apply fft and mv-op
    fJ = fftn(JIn0(:,:,:,1),[LfG, MfG, NfG]);
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn0(:,:,:,2),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn0(:,:,:,3),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);
    
    % sum contribution and return local indexes
    t = dx3 * s - JOut(idx);
    % t = A*s; % apply the operation on s
    % ---------------------------------------------------------------------
    
    omega = (t'*s)/(t'*t); % compute omega
    x = x + omega*s; % update solution
    r = s - omega*t; % update residual
    currentii = 2*(ii+1)+1;
    resvec(currentii) = norm(r)/bnorm; % store half iteration residual
    
    % check if better solution that current one
    if resvec(currentii) < resvec(bestii)
        
        bestii = currentii;
        
        % if better solution check half convergence
        if ( resvec(bestii) <= tol)
            flag = 0; % converged
            break;
        end
        
        % if not, store as best solution so far
        bestsol = x;
        stag = 0;

        
    else
        
        % check for stagnation or bad conditioning of variables
        if abs(omega)*resvec(currentii-1) < 1e-10
            if (stag > 5)
                flag = 3; % stagnated
                break
            else
                stag = stag + 1;
            end
        else
            if isfinite(omega)
                stag = 0;
            else
                flag = 4; % numerical errors
                break
            end
        end
        
    end
        
    % continue: update data for next iteration
    rhonew = ro'*r; % new rho   
    beta = rhonew*alpha/(rho*omega); % beta
    p = r + beta*(p - omega*v); % update p for next iteration
    rho = rhonew; % assign rho
    
end

      
% check final values and return data
resvec = resvec(1:currentii);
resvec = gather(resvec);

switch flag
    
    case 0 % converged, x is final result
        x = gather(x);
        relres = resvec(currentii);
        iter = [currentii currentii];
        
    case 3 % stagnated, not converged
        x = gather(bestsol);
        relres = resvec(bestii);
        iter = [currentii bestii];

    case 4 % numerical errors
        x = gather(bestsol);
        relres = resvec(bestii);
        iter = [currentii bestii];
        
    otherwise % end of iterations
        if (bestii < currentii) % best solution is not the final one
            x = gather(bestsol);
            relres = resvec(bestii);
            iter = [currentii bestii];
            if (relres <= tol)
                flag = 0; % convergence
            else
                flag = 1; % no convergence
            end
        else % best solution is the final one (may differ from bestsol)
            x = gather(x);
            relres = resvec(currentii);
            iter = [currentii currentii];
            if (relres <= tol)
                flag = 0; % convergence
            else
                flag = 1; % max iterations, no convergence
            end
        end
        
end
