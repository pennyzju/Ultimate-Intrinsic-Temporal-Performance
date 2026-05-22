function [Eout,Hout] = dda_HEradiateGPU(Jin, fG, fK, k0, L, M ,N)
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Amit's code for free space radiation
%
% _________________________________________________________________________
% _________________________________________________________________________



% reshape and initialize the arrays and transfer from local to global
Eout = gpuArray(Jin); % copy vector in Eout
Eout = reshape(Eout,L,M,N,3);
Hout = gpuArray.zeros(L, M, N, 3);

% apply fft and mv-op for each component and add to the corresponding one

% ------ 1st component
fJ = fftn(Eout(:,:,:,1),[2*L, 2*M, 2*N]);
% contribution for E field
Eout1 = fG(:,:,:,1,1) .* fJ;
Eout2 = fG(:,:,:,2,1) .* fJ;
Eout3 = fG(:,:,:,3,1) .* fJ;
% contribution for H field
% Hout1 is zero
Hout2 = -fK(:,:,:,3) .* fJ; % second component of Hout = -z Jin1
Hout3 = fK(:,:,:,2) .* fJ;  % third component of Hout = y Jin1

% ------ 2nd component
fJ = fftn(Eout(:,:,:,2),[2*L, 2*M, 2*N]);
% contribution for E field
Eout1 = Eout1 + fG(:,:,:,1,2) .* fJ;
Eout2 = Eout2 + fG(:,:,:,2,2) .* fJ;
Eout3 = Eout3 + fG(:,:,:,3,2) .* fJ;
% contribution for H field
Hout1 = fK(:,:,:,3) .* fJ; % first component of Hout = +z Jin2
% Hout2 is zero
Hout3 = Hout3 - fK(:,:,:,1) .* fJ;  % third component of Hout = -x Jin2

% ------ 3rd component
fJ = fftn(Eout(:,:,:,3),[2*L, 2*M, 2*N]);
% contribution for E field
Eout1 = Eout1 + fG(:,:,:,1,3) .* fJ;
Eout2 = Eout2 + fG(:,:,:,2,3) .* fJ;
Eout3 = Eout3 + fG(:,:,:,3,3) .* fJ;
% contribution for H field    
Hout1 = Hout1 - fK(:,:,:,2) .* fJ; % first component of Hout = -y Jin3
Hout2 = Hout2 + fK(:,:,:,1) .* fJ;  % third component of Hout = +x Jin3
% Hout3 is zero

% apply ifft and add contribution of the incident
eta =  3.767303134617706e+002; % Free-space impedance
mult = (3*1j*k0/eta);
Eout1 = ifftn(Eout1);
Eout(:,:,:,1) = Eout1(L:2*L-1, M:2*M-1, N:2*N-1) - Eout(:,:,:,1)/mult;
Eout2 = ifftn(Eout2);
Eout(:,:,:,2) = Eout2(L:2*L-1, M:2*M-1, N:2*N-1) - Eout(:,:,:,2)/mult;
Eout3 = ifftn(Eout3);
Eout(:,:,:,3) = Eout3(L:2*L-1, M:2*M-1, N:2*N-1) - Eout(:,:,:,3)/mult;
       
% apply ifft and add contribution of the incident
Hout1 = ifftn(Hout1);
Hout(:,:,:,1) = Hout1(L:2*L-1, M:2*M-1, N:2*N-1);
Hout2 = ifftn(Hout2);
Hout(:,:,:,2) = Hout2(L:2*L-1, M:2*M-1, N:2*N-1);
Hout3 = ifftn(Hout3);
Hout(:,:,:,3) = Hout3(L:2*L-1, M:2*M-1, N:2*N-1);
   
% get data back from GPU
Eout = gather(Eout);
Hout = gather(Hout);


