function [Hout] = dda_Hradiate(Jin, fK, L, M ,N)
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Amit's code for free space radiation
%
% _________________________________________________________________________
% _________________________________________________________________________

% reshape and initialize the arrays and transfer from local to global
Jin0 = zeros(L, M, N, 3);
Jin0(:) = Jin;
Hout = zeros(L, M, N, 3);

% [Lk,Mk,Nk,~] = size(fK);


% apply fft and mv-op for each component and add to the corresponding one

%first component of Jin
fJ = fftn(Jin0(:,:,:,1),[2*L, 2*M, 2*N]);
% Jout1 is zero
Jout2 = -fK(:,:,:,3) .* fJ; % second component of Jout = -z Jin1
Jout3 = fK(:,:,:,2) .* fJ;  % third component of Jout = y Jin1

%second component of Jin
fJ = fftn(Jin0(:,:,:,2),[2*L, 2*M, 2*N]);
Jout1 = fK(:,:,:,3) .* fJ; % first component of Jout = +z Jin2
% Jout2 is zero
Jout3 = Jout3 - fK(:,:,:,1) .* fJ;  % third component of Jout = -x Jin2

%third component of Jin
fJ = fftn(Jin0(:,:,:,3),[2*L, 2*M, 2*N]);
Jout1 = Jout1 - fK(:,:,:,2) .* fJ; % first component of Jout = -y Jin3
Jout2 = Jout2 + fK(:,:,:,1) .* fJ;  % third component of Jout = +x Jin3
% Jout3 is zero

    
% apply ifft and add contribution of the incident
Jout1 = ifftn(Jout1);
Hout(:,:,:,1) = Jout1(L:2*L-1, M:2*M-1, N:2*N-1);
Jout2 = ifftn(Jout2);
Hout(:,:,:,2) = Jout2(L:2*L-1, M:2*M-1, N:2*N-1);
Jout3 = ifftn(Jout3);
Hout(:,:,:,3) = Jout3(L:2*L-1, M:2*M-1, N:2*N-1);
   


