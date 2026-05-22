close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
lambda  = 2*pi;
ko = 2*pi/lambda;
f = co/lambda;
omega = 2 * pi * f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ko_L = 1.0; % x = k*a
%
Length = ko_L / ko;
%
epsilon = 25.0 - 1i*30.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nX = 51 ; % Number of cells along each direction
% Define grid
d = 1.15*(Length/2);
dx = 2*d/nX;
dy = dx;
dz = dx;
%
x = -d + ((1:nX)-0.5)*dx;
%
r = grid3d(x);
% Electric properties
cube = @(r) ( abs(r(:,:,:,1)) < (Length/2) ) & ( abs(r(:,:,:,2)) < (Length/2) ) & (abs( r(:,:,:,3) ) < (Length/2));
fEpsilon = @(r)cube(r)*(epsilon-1)+1;
e_r = fEpsilon(r);
%
name_er = sprintf('./DATA/er_CUBE_Vox%d_e%d_%d.mat',nX,real(epsilon),abs(imag(epsilon)));
save(name_er, 'e_r')