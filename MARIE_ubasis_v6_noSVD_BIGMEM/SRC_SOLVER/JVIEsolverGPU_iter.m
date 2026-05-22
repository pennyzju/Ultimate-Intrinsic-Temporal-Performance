% GPU adapted Preconditioned GMRES internal iteration routine for JVIE solver
%
%   Version by J. Fernandez Villena - jvillena@mit.edu
%              Computational Prototyping Group, RLE at MIT
%
% Solves A z = r for z, then returns x + z
%
% INPUT:  x      current solution vector
%         r      N-by-1 preconditioned residual vector
%         m      number of GMRES iterations to perform
%         tol    specifies the tolerance of the method
%         fG        Circulant matrix for function application
%         tau       inhomogeneous properties at given frequency
%         L         number of discretization elements in x direction
%         M         number of discretization elements in y direction
%         N         number of discretization elements in z direction
%         dx3       volume of discretized voxel
%         idx       indexes of the elements relevant for the solution
%
% OUTPUT: x      updated solution vector
%         r      preconditioned residual vector
%         k      number of GMRES iterations actually performed
%         resvec vector containing norm of residual at each iteration of GMRES
%
%


function [x,r,k,resvec] = JVIEsolverGPU_iter(x,r,m,tol,fG,tau,L,M,N,dx3,idx)

% compute norm or current residue
rnorm = norm(r);

% Initialize V
V = gpuArray.zeros(size(r,1),m+1);
V(:,1) = r / rnorm;

% Preallocate H and resvec
H = gpuArray.zeros(m+1,m);
resvec = gpuArray.zeros(m,1);

% initialize function data
JIn0 = gpuArray.zeros(L, M, N, 3);
JOut = gpuArray.zeros(L, M, N, 3);

% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);


for k = 1:m
    
    w = V(:,k);
    
    % -------------------------------------------------------------------------
    %   Apply the function A*w
    %
    JIn0(idx) = w;

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
    w = dx3 * w - JOut(idx);
   
    % -------------------------------------------------------------------------
    %   Create next column of V and H
    %
    
    H(1:k,k) = V(:,1:k)' * w;
    w = w - V(:,1:k)*H(1:k,k);
    
    H(k+1,k) = norm(w);
    V(:,k+1) = w / H(k+1,k);
    
    % Initialize right hand side of least-squares system
    rhs = gpuArray.zeros(k+1,1);
    rhs(1) = norm(r);
    
    % Solve least squares system; Calculate residual norm
    Hc = H(1:k+1,1:k);
    y = Hc \ rhs;
    res = rhs - Hc * y;
    resvec(k) = norm(res);
    
    % check for early convergence
    if resvec(k) < tol
        
        break;
        
    end
    
end

% Calculate solution and residual
resvec = resvec(1:k);
% apply update x = x + V(:,1:k) * y;
x = x + V(:,1:k)*y;
r = V(:,1:k+1)*res;
       