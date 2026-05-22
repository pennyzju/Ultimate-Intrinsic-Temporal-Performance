function [Hout] = H_field_Kop_compGPU(J, fK, dV, Hinc)
%
% Function to compute the H field
%   Computational Prototyping Group, RLE at MIT
%

% fft dimensions
[LfK, MfK, NfK, ~] = size(fK);

% vector dimensions
[L, M, N, ~] = size(Hinc);
J = reshape(J, L, M, N, 3);

% apply fft and mv-op
J1 = gpuArray( J(:,:,:,1) );
fJ = gather( fftn( J1,[LfK, MfK, NfK]) ); % first component Vin
clear J1;
Jout3 = fK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x
Jout2 = -fK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
% Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x

J2 = gpuArray( J(:,:,:,2) );
fJ = gather( fftn(J2,[LfK, MfK, NfK]) ); % second component Vin
clear J2;
Jout1 = fK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
% Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
Jout3 = Jout3 - fK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y

J3 = gpuArray( J(:,:,:,3) );
fJ = gather( fftn(J3,[LfK, MfK, NfK]) ); % third component Vin
clear J3;
% Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z
Jout2 = Jout2 + fK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
Jout1 = Jout1 - fK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z


% apply ifft ( Hout = (Hinc+Jout)/dV )
% Hout = Hinc ./dV;
Hout = Hinc;   % *** BUG FIX BY THANOS ***

Jout1 = gather( ifftn( gpuArray(Jout1) ) );
Hout(:,:,:,1) = Hout(:,:,:,1) + Jout1(1:L,1:M,1:N) ./ dV;

Jout2 = gather( ifftn( gpuArray(Jout2) ) );
Hout(:,:,:,2) = Hout(:,:,:,2) + Jout2(1:L,1:M,1:N) ./ dV;

Jout3 = gather( ifftn( gpuArray(Jout3) ) );
Hout(:,:,:,3) = Hout(:,:,:,3) + Jout3(1:L,1:M,1:N) ./ dV;




