function [E] = E_field_Nop_compGPU_bastien(J, fG,  dV, omega, eo, Einc)



[L, M, N, ~] = size(Einc);

J = reshape(J, L, M, N, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);


% Compute FFT(J) and mv-op
J1 = gpuArray( J(:,:,:,1) );
fJ = gather( fftn(J1, [LfG, MfG, NfG]) );
clear J1;
Jout1 = fG(:,:,:,1) .* fJ;
Jout2 = fG(:,:,:,2) .* fJ;
Jout3 = fG(:,:,:,3) .* fJ;

J2 = gpuArray( J(:,:,:,2) );
fJ = gather( fftn(J2, [LfG, MfG, NfG]) );
clear J2;
Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
Jout3 = Jout3 + fG(:,:,:,5) .* fJ;

J3 = gpuArray( J(:,:,:,3) );
fJ = gather( fftn(J3, [LfG, MfG, NfG]) );
clear J3;
Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
Jout3 = Jout3 + fG(:,:,:,6) .* fJ;

% apply ifft
fJ = gather( ifftn( gpuArray(Jout1) ) );
E(:,:,:,1) = fJ(1:L, 1:M, 1:N);
fJ = gather( ifftn( gpuArray(Jout2) ) );
E(:,:,:,2) = fJ(1:L, 1:M, 1:N);
fJ = gather( ifftn( gpuArray(Jout3) ) );
E(:,:,:,3) = fJ(1:L, 1:M, 1:N);

% assembly field
E = E - dV * J;
E = (1/dV) * E;

% E = (1/dV) * Einc + E./(1i*omega*eo); % *** BUG FIX BY THANOS ***
E = Einc + E./(1i*omega*eo); 

















