% GPU adapted Dedicated Solver for JVIE formulation
% using restarted GMRES 
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
%         restart   inner iterations before restart
%         tol       relative tolerance for the residue
%         maxit     maximum number of outer iterations (total is restart*maxit)
%         x0        initial guess
%
% OUTPUT: x         solution vector
%         flag      0 if converged, 1 if not
%         relres    final relative residue: norm(b - A*x)/norm(b) 
%         iter      vector with [current internal iterations, current external iterations]
%         resvec    vector containing norm of residual at each iteration of GMRES
%         
%
%
% SUBFUNCTIONS: fAgpu_d2, JVIEsolverGPU_iter.m
%
%

function [x,flag,relres,iter,resvec] = JVIEsolverGPU(fG, tau, L, M, N, dx3, idx, b, restart, tol, maxit, x0)

% Initialize variables.
if(nargin < 8 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 9 || isempty(restart))
   restart = 50;
end
if(nargin < 10 || isempty(tol))
   tol = 1e-3;
end
if(nargin < 11 || isempty(maxit))
   maxit = 10;
end
if(nargin < 12 || isempty(x0))
    x0 = gpuArray.zeros(size(b));
end


% -------------------------------------------------------------------------
% settings

% send data to GPU
fG = gpuArray(fG);
tau = gpuArray(tau);
idx = gpuArray(idx);
b = gpuArray(b);
x = gpuArray(x0);


% Calculate initial preconditioned residual.
% r = b - A*x0;
r = b - fAgpu_d2(x0,fG,tau,L,M,N,dx3,idx);

% Calculate rhs norm
bnorm = norm(b);

% initialize the residue vector and other variables
resvec = gpuArray.zeros(restart*maxit,1);
it = 1;
outit = 0;
resvec(it) = norm(r)/bnorm;


% -------------------------------------------------------------------------
% Restarted loop until convergence or maximum reached
while(resvec(it) > tol) && (outit < maxit)
  
         
    % Call the gmres to perform restart iterations ------------------------
    [x,r,p,resvec_inner] = JVIEsolverGPU_iter(x,r,restart,tol*bnorm,fG,tau,L,M,N,dx3,idx);
           
    % store vector with residues
    resvec(it+1:it+p) = resvec_inner/bnorm;
    it = it + p;
    
    % %     % Calculate relative residual.
    % %     relres = norm(b -gpufA_d2(x,fG,tau,L,M,N,dx3))/bnorm;
    % %     fprintf(1,'restart %d, it %d, ||r|| = %g ( resvec %g, tol %g, reltol %g)\n' ,outit, it, relres, resvec(it), tol, tol*bnorm);
    
    % increase external loop counter
    outit = outit+1;

   
end

% Calculate residual, and iteration count
resvec = resvec(1:it);
iter = [p outit];

% Calculate relative residual.
relres = norm(b - fAgpu_d2(x,fG,tau,L,M,N,dx3,idx))/bnorm;
x = gather(x);
resvec = gather(resvec);
if (relres <= tol ) % converged
    flag = 0;
else
    flag = 1;
end

