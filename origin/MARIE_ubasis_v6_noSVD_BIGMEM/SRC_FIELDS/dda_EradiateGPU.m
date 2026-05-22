function [Eout] = dda_EradiateGPU(Jin, fG, k0, L, M ,N)
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Amit's code for free space radiation
%
% _________________________________________________________________________
% _________________________________________________________________________



% reshape and initialize the arrays and transfer from local to global
Jin0 = gpuArray(Jin);
Jin0 = reshape(Jin0,L,M,N,3);
Eout = gpuArray.zeros(L, M, N, 3);

% apply fft and mv-op for each component and add to the corresponding one
fJ = fftn(Jin0(:,:,:,1),[2*L, 2*M, 2*N]);
Jout1 = fG(:,:,:,1,1) .* fJ;
Jout2 = fG(:,:,:,2,1) .* fJ;
Jout3 = fG(:,:,:,3,1) .* fJ;

fJ = fftn(Jin0(:,:,:,2),[2*L, 2*M, 2*N]);
Jout1 = Jout1 + fG(:,:,:,1,2) .* fJ;
Jout2 = Jout2 + fG(:,:,:,2,2) .* fJ;
Jout3 = Jout3 + fG(:,:,:,3,2) .* fJ;

fJ = fftn(Jin0(:,:,:,3),[2*L, 2*M, 2*N]);
Jout1 = Jout1 + fG(:,:,:,1,3) .* fJ;
Jout2 = Jout2 + fG(:,:,:,2,3) .* fJ;
Jout3 = Jout3 + fG(:,:,:,3,3) .* fJ;
    
% apply ifft and add contribution of the incident
eta =  3.767303134617706e+002; % Free-space impedance
mult = (3*1j*k0/eta);
Jout1 = ifftn(Jout1);
Eout(:,:,:,1) = Jout1(L:2*L-1, M:2*M-1, N:2*N-1) - Jin0(:,:,:,1)/mult;
Jout2 = ifftn(Jout2);
Eout(:,:,:,2) = Jout2(L:2*L-1, M:2*M-1, N:2*N-1) - Jin0(:,:,:,2)/mult;
Jout3 = ifftn(Jout3);
Eout(:,:,:,3) = Jout3(L:2*L-1, M:2*M-1, N:2*N-1) - Jin0(:,:,:,3)/mult;
   
% get data back from GPU
Eout = gather(Eout);


