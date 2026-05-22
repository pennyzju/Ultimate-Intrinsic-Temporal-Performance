function [fG] = FFT_Gmn(Gp_mn)
% size of Gp_mn
% (2*L-1,2*M-1,2*N-1)
[Lfft, Mfft, Nfft, ~] = size(Gp_mn);


% 6 block case
fG = zeros(Lfft, Mfft, Nfft, 6);
% 3D-FFT Gp_mn

fG(:,:,:,1) = fftn(Gp_mn(:,:,:,1));
fG(:,:,:,2) = fftn(Gp_mn(:,:,:,2));
fG(:,:,:,3) = fftn(Gp_mn(:,:,:,3));
fG(:,:,:,4) = fftn(Gp_mn(:,:,:,4));
fG(:,:,:,5) = fftn(Gp_mn(:,:,:,5));
fG(:,:,:,6) = fftn(Gp_mn(:,:,:,6));

