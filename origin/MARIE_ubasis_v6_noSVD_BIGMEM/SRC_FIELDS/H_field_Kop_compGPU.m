function [Hout] = H_field_Kop_compGPU(J, fK, dV, Hinc)
%
% Function to compute the H field
%   Computational Prototyping Group, RLE at MIT
%

% fft dimensions
[LfK, MfK, NfK, ~] = size(fK);

% vector dimensions
[L, M, N, ~] = size(Hinc);
Hinc = gpuArray(Hinc);
J = gpuArray(J);
J = reshape(J, L, M, N, 3);

% allocate space
Hout = gpuArray.zeros(L, M, N, 3);


% apply fft and mv-op
fJ = fftn(J(:,:,:,1),[LfK, MfK, NfK]); % first component Vin
Jout3 = fK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x
Jout2 = -fK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
% Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x

fJ = fftn(J(:,:,:,2),[LfK, MfK, NfK]); % second component Vin
Jout1 = fK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
% Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
Jout3 = Jout3 - fK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y

fJ = fftn(J(:,:,:,3),[LfK, MfK, NfK]); % third component Vin
% Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z
Jout2 = Jout2 + fK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
Jout1 = Jout1 - fK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z


% apply ifft
Jout1 = ifftn(Jout1);
Hout(:,:,:,1) = Jout1(1:L,1:M,1:N);
Jout2 = ifftn(Jout2);
Hout(:,:,:,2) = Jout2(1:L,1:M,1:N);
Jout3 = ifftn(Jout3);
Hout(:,:,:,3) = Jout3(1:L,1:M,1:N);



% assembly field
% Hout = (Hinc + Hout)./dV;     % *** BUG FIX BY THANOS ***
Hout = Hinc  +  Hout./dV;

Hout = gather(Hout); 
