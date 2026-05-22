function [JOut] = mv_op_fft_Nop_2_semigpu(JIn0, fG, e_r, dV, transp_flag, idx)
%
% Function that applies the JVIE II core operation
%   Computational Prototyping Group, RLE at MIT
%

% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);

[L, M, N] = size(e_r);

% allocate space
Jgpu = gpuArray.zeros(L, M, N, 3);
fGgpu = gpuArray.zeros(LfG, MfG, NfG);

% translate from local to global coordinates
Jgpu(idx) = gpuArray(JIn0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(transp_flag,'transp')       % y = A'*x
    
    tau = conj((e_r - 1.0) ./ e_r);
    tau = gpuArray(tau);
    
    % apply fft and mv-op
    fJ = fftn(conj( tau .* Jgpu(:,:,:,1)) , [LfG, MfG, NfG]);    
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,1));
    Jout1 = fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,3));
    Jout3 = fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,2));
    Jout2 = fGgpu .* fJ;
    
    fJ = fftn(conj( tau .* Jgpu(:,:,:,2)) , [LfG, MfG, NfG]);
    % fGgpu(:,:,:) = gpuArray(fG(:,:,:,2));
    Jout1 = Jout1 + fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,4));
    Jout2 = Jout2 + fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,5));
    Jout3 = Jout3 + fGgpu .* fJ;
    
    fJ = fftn(conj( tau .* Jgpu(:,:,:,3)) , [LfG, MfG, NfG]);
    % fGgpu(:,:,:) = gpuArray(fG(:,:,:,5));
    Jout2 = Jout2 + fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,6));
    Jout3 = Jout3 + fGgpu .* fJ;
    fGgpu(:,:,:) = gpuArray(fG(:,:,:,3));
    Jout1 = Jout1 + fGgpu .* fJ;
    
    clear fGgpu; clear fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    Jgpu(:,:,:,1) = conj( Jout1(1:L,1:M,1:N));
    Jout2 = ifftn(Jout2);
    Jgpu(:,:,:,2) = conj( Jout2(1:L,1:M,1:N));
    Jout3 = ifftn(Jout3);
    Jgpu(:,:,:,3) = conj( Jout3(1:L,1:M,1:N));
    
elseif strcmp(transp_flag,'notransp') % y = A*x
    
    tau = (e_r - 1.0) ./ e_r ;
   
    % apply fft and mv-op
    fJ = fftn(Jgpu(:,:,:,1),[LfG, MfG, NfG]);
    fGgpu(:,:,:) = fG(:,:,:,1);
    Jout1 = fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,3);
    Jout3 = fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,2);
    Jout2 = fGgpu .* fJ;
    
    fJ = fftn(Jgpu(:,:,:,2),[LfG, MfG, NfG]);
    % fGgpu(:,:,:) = fG(:,:,:,2);
    Jout1 = Jout1 + fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,4);
    Jout2 = Jout2 + fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,5);
    Jout3 = Jout3 + fGgpu .* fJ;
    
    fJ = fftn(Jgpu(:,:,:,3),[LfG, MfG, NfG]);
    % fGgpu(:,:,:) = fG(:,:,:,5);
    Jout2 = Jout2 + fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,6);
    Jout3 = Jout3 + fGgpu .* fJ;
    fGgpu(:,:,:) = fG(:,:,:,3);
    Jout1 = Jout1 + fGgpu .* fJ;

    clear fGgpu; clear fJ;
    tau = gpuArray(tau);
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    Jgpu(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    Jgpu(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    Jgpu(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);

end

JOut = gather(Jgpu);

% return local coordinates
JOut = dV * JIn0 - JOut(idx);

clear Jgpu; clear tau; clear Jout1; clear Jout2; clear Jout3;