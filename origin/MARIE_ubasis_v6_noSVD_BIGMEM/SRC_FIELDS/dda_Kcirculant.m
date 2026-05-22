function [fK] = dda_Kcirculant(r, dx, ko)
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Amit's code for generating the circulant
%
% _________________________________________________________________________
% _________________________________________________________________________



% get grid size
[L, M, N, ~] = size(r);

% compute second grid
x = ((1:2*L-1)-L)*dx; 
y = ((1:2*M-1)-M)*dx;
z = ((1:2*N-1)-N)*dx;
r2 = grid3d(x, y, z);
clear x; clear y; clear z;

R = dda_abs3d(r2);
% unitR = dda_scalar_mult(1./R, r2);
R3 = R.^3;

fK = zeros(2*L, 2*M, 2*N, 3);

% Calc dyadic Green's function
mult = (1 + 1j*ko*R).*exp(-1j*ko*R)./(4*pi*R3);

Kx = mult.*r2(:,:,:,1);
Kx(L,M,N) = 0;
Ky = mult.*r2(:,:,:,2);
Ky(L,M,N) = 0;
Kz = mult.*r2(:,:,:,3);
Kz(L,M,N) = 0;

fK(:,:,:,1) = fftn(Kx, [2*L, 2*M, 2*N]);
fK(:,:,:,2) = fftn(Ky, [2*L, 2*M, 2*N]);
fK(:,:,:,3) = fftn(Kz, [2*L, 2*M, 2*N]);


