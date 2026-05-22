% Function that applies the JVIE core operation
%   Computational Prototyping Group, RLE at MIT
%
% INPUT:  JIn0      excitation vector
%         fG        Circulant matrix for function application
%         tau       inhomogeneous properties at given frequency
%         L         number of discretization elements in x direction
%         M         number of discretization elements in y direction
%         N         number of discretization elements in z direction
%         dx3       volume of discretized voxel
%         idx       indexes of the elements relevant for the solution
%
% OUTPUT: JOut      solution vector
%
%

function [JOut] = fA_d2(JIn0,fG,tau,L,M,N,dx3,idx)


% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);

% initialize function data
JIn = zeros(L, M, N, 3);
JOut = zeros(L, M, N, 3);

% transform local data to global
JIn(idx) = JIn0;

% apply fft and mv-op
fJ = fftn(JIn(:,:,:,1),[LfG, MfG, NfG]);
Jout1 = fG(:,:,:,1) .* fJ;
Jout2 = fG(:,:,:,2) .* fJ;
Jout3 = fG(:,:,:,3) .* fJ;

fJ = fftn(JIn(:,:,:,2),[LfG, MfG, NfG]);
Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
Jout3 = Jout3 + fG(:,:,:,5) .* fJ;

fJ = fftn(JIn(:,:,:,3),[LfG, MfG, NfG]);
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

%return local coordinates
JOut = dx3 * JIn0 - JOut(idx);