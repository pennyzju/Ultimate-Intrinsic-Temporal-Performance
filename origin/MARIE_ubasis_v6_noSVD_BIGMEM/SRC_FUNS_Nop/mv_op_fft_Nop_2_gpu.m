function [JOut] = mv_op_fft_Nop_2_gpu(JIn0, fG, e_r, dV, transp_flag, idx)
%
% GPU accelerated function that applies the JVIE II core operation
%   Computational Prototyping Group, RLE at MIT
%
% fG and e_r must be in GPU mem


% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);

[L, M, N] = size(e_r);

% allocate space
JOut = gpuArray.zeros(L, M, N, 3);
JIn = gpuArray.zeros(L, M, N, 3);

% send to gpu
JIn(idx) = gpuArray(JIn0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(transp_flag,'transp')       % y = A'*x
    
    tau = conj((e_r - 1.0) ./ e_r);
    
    % apply fft and mv-op
    fJ = fftn(conj( tau .* JIn(:,:,:,1)) , [LfG, MfG, NfG]);
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = fftn(conj( tau .* JIn(:,:,:,2)) , [LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = fftn(conj( tau .* JIn(:,:,:,3)) , [LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = conj( Jout1(1:L,1:M,1:N));
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = conj( Jout2(1:L,1:M,1:N));
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = conj( Jout3(1:L,1:M,1:N));
    
elseif strcmp(transp_flag,'notransp') % y = A*x
    
    tau = (e_r - 1.0) ./ e_r;

    % apply fft and mv-op
    fJ = fftn(JIn(:,:,:,1),[LfG, MfG, NfG]);
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn(:,:,:,2),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn(:,:,:,3),[LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);    
    JOut(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);
    
end

% apply op and return local coordinates
JOut = dV * JIn(idx) - JOut(idx);

% get from GPU
JOut = gather(JOut);

% clear gpu data
clear JIn; clear tau; clear Jout1; clear Jout2; clear Jout3; clear fJ;