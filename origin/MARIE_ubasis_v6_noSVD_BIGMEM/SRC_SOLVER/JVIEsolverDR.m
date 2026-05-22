
% Dedicated Solver for JVIE formulation, using GMRES DR
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
%         ritz      number of Ritz singular values to use for deflation
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
% SUBFUNCTIONS: fA_d2.m, JVIEsolverDR_iter.m
%
%


function [x,flag,relres,iter,resvec] = JVIEsolverDR(fG, tau, L, M, N, dx3, idx, b, restart, tol, maxit, ritz, x0)

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
if(nargin < 12 || isempty(ritz))
    ritz = min(10,round(restart/2));
else
    ritz = min(ritz,restart-1);
end
if(nargin < 13 || isempty(x0))
    x0 = zeros(size(b));
end



% -------------------------------------------------------------------------
% settings

% Calculate initial preconditioned residual.
% r = b - A*x0;
r = b - fA_d2(x0,fG,tau,L,M,N,dx3,idx);

% Calculate rhs norm
bnorm = norm(b);

% initialize the residue vector and other variables
resvec = zeros(restart*maxit,1);
it = 1;
outit = 0;
resvec(it) = norm(r)/bnorm;
C = [];


% -------------------------------------------------------------------------
% first call to an initial gmres to perform restart iterations

% Call the gmres to perform restart iterations
[V,H,~,p,resvec_inner] = JVIEsolverDR_iter(r,restart,C,tol*bnorm,fG,tau,L,M,N,dx3,idx);

% store vector with residues
resvec(it+1:it+p) = resvec_inner/bnorm;
it = it + p;

% obtain update on solution and residual
y = H\(V'*r);
x = V(:,1:p)*y;
r = r - V*(H*y);


% % % Calculate relative residual.
% % relres = norm(b - gpu_mv_op_fft_d2(x0+x,fG,tau,L,M,N,dx3))/norm(b);
% % fprintf(1,'restart %d, it %d, ||r|| = %g ( resvec %g, tol %g, reltol %g)\n' ,outit, it, relres, resvec(it), tol, tol*bnorm);


% Generate the Ritz vectors for deflation ----------------------------
if (resvec(it) > tol)

    % Find the ritz smallest harmonic Ritz vectors.
    if ritz > p-1
        P = getHarmVecs1(p,p-1,H);
    else
        P = getHarmVecs1(p,ritz,H);
    end
    
    % Form the subspace to recycle: U
    U = V(:,1:p)*P;
    
    % Form orthonormalized C
    [C,R] = qr(H*P,0);
    C = V*C;
    
    % adjust U accordingly so that C = A*U
    U = U/R;
    
    % clean stuff
    clear V; clear H; clear P;
    
end


% -------------------------------------------------------------------------
% Restarted loop until convergence or maximum reached
while(resvec(it) > tol) && (outit < maxit)
  
    % increase external loop counter
    outit = outit+1;
        
    % Call the gmres to perform restart iterations ------------------------
    [V,H,B,p,resvec_inner] = JVIEsolverDR_iter(r,restart,C,tol*bnorm,fG,tau,L,M,N,dx3,idx);
        
    % store vector with residues
    resvec(it+1:it+p) = resvec_inner/bnorm;
    it = it + p;
           
   % Rescale U; 
   % Store inverses of the norms of columns of U in diagonal matrix D
   d = zeros(ritz,1);
   for i = 1:ritz
      d(i) = norm(U(:,i));
   end
   D = diag(1./d);
   U = U*D;

   % Form large H
   H2 = [D, B; zeros(size(H,1),size(D,2)), H];

   % obtain update on solution and residual
   y = H2\([C V]'*r);
   x = x + [U V(:,1:p)]*y;
   r = r - [C V]*(H2*y);
   
   % %    % Calculate relative residual.
   % %    relres = norm(b - gpu_mv_op_fft_d2(x0+x,fG,tau,L,M,N,dx3))/norm(b);
   % %    fprintf(1,'restart %d, it %d, ||r|| = %g ( resvec %g, tol %g, reltol %g)\n' ,outit, it, relres, resvec(it), tol, tol*bnorm);

   
   if (resvec(it) > tol)
       % Calculate Harmonic Ritz vectors.
       P = getHarmVecs2(p+ritz,ritz,H2,V,U,C);
       
       % Form new U and C.
       U = [U V(:,1:p)]*P;
       
       % Form orthonormalized C and adjust U accordingly so that C = A*U
       [Q,R] = qr(H2*P,0);
       C = [C V] * Q;
       U = U/R;
   end
   
   
end

% Calculate final solution,
x = x0 + x;

% Calculate residual, and iteration count
resvec = resvec(1:it);
iter = [p outit];

% Calculate relative residual.
relres = norm(b - fA_d2(x,fG,tau,L,M,N,dx3,idx))/norm(b);
if (relres <= tol ) % converged
    flag = 0;
else
    flag = 1;
end

