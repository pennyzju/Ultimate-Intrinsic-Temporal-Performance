function [Vout] = mv_op_fft_Kop(Vin, fK)
%
% Function that applies the JVIE core Kop operation
%   Computational Prototyping Group, RLE at MIT
%

% fft dimensions
[LfK, MfK, NfK, ~] = size(fK);

% vector dimensions
[L, M, N, ~] = size(Vin);

% allocate space
Vout = zeros(L, M, N, 3);


% apply fft and mv-op
fJ = fftn(Vin(:,:,:,1),[LfK, MfK, NfK]); % first component Vin
Jout3 = fK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x
Jout2 = -fK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
% Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x

fJ = fftn(Vin(:,:,:,2),[LfK, MfK, NfK]); % second component Vin
Jout1 = fK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
% Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
Jout3 = Jout3 - fK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y

fJ = fftn(Vin(:,:,:,3),[LfK, MfK, NfK]); % third component Vin
% Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z
Jout2 = Jout2 + fK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
Jout1 = Jout1 - fK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z


% apply ifft
Jout1 = ifftn(Jout1);
Vout(:,:,:,1) = Jout1(1:L,1:M,1:N);
Jout2 = ifftn(Jout2);
Vout(:,:,:,2) = Jout2(1:L,1:M,1:N);
Jout3 = ifftn(Jout3);
Vout(:,:,:,3) = Jout3(1:L,1:M,1:N);
   

