function [fG] = FFT_CK(Kop_mn)
% size of Gp_mn
% (2*L-1,2*M-1,2*N-1)
[Lfft, Mfft, Nfft, ~] = size(Kop_mn);
%
fG = zeros(Lfft, Mfft, Nfft, 3);
% 3D-FFT Gp_mn
for p=1:3
    fG(:,:,:,p) = fftn(Kop_mn(:,:,:,p));
end