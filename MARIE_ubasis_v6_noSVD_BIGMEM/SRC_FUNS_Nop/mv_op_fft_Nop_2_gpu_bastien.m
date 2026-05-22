function [JOut] = mv_op_fft_Nop_2_gpu_bastien(JIn0, fG, e_r, dV, transp_flag, idx)
%
% GPU accelerated function that applies the JVIE II core operation
%   Computational Prototyping Group, RLE at MIT
%


% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);

[L, M, N] = size(e_r);

% allocate space
JOut = zeros(L, M, N, 3);

JIn  = zeros(L, M, N, 3);
JIn(idx) = JIn0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(transp_flag,'transp')       % y = A'*x
    
    tau = conj((e_r - 1.0) ./ e_r);
    
    % apply fft and mv-op
    tmp = gpuArray( tau .* JIn(:,:,:,1) );
    fJ = gather( fftn(conj(tmp) , [LfG, MfG, NfG]) );
    clear tmp;
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    tmp = gpuArray( tau .* JIn(:,:,:,2) );
    fJ = fftn(conj(tmp) , [LfG, MfG, NfG]);
    clear tmp;
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    tmp = gpuArray( tau .* JIn(:,:,:,3) );
    fJ = fftn(conj(tmp) , [LfG, MfG, NfG]);
    clear tmp;
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = gather( ifftn( gpuArray(Jout1 ) ) );
    JOut(:,:,:,1) = conj( Jout1(1:L,1:M,1:N));
    
    Jout2 = gather( ifftn( gpuArray(Jout2) ) );
    JOut(:,:,:,2) = conj( Jout2(1:L,1:M,1:N));
    
    Jout3 = gather( ifftn( gpuArray(Jout3) ) );
    JOut(:,:,:,3) = conj( Jout3(1:L,1:M,1:N));
    
elseif strcmp(transp_flag,'notransp') % y = A*x
    
    tau = (e_r - 1.0) ./ e_r;

    % apply fft and mv-op
    fJ = gather( fftn( gpuArray(JIn(:,:,:,1)),[LfG, MfG, NfG]) );
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = gather( fftn( gpuArray(JIn(:,:,:,2)),[LfG, MfG, NfG]) );
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = gather( fftn( gpuArray(JIn(:,:,:,3)),[LfG, MfG, NfG]) );
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = gather( ifftn( gpuArray(Jout1) ) );
    JOut(:,:,:,1) = tau .* Jout1(1:L,1:M,1:N);
    
    Jout2 = gather( ifftn( gpuArray(Jout2) ) );
    JOut(:,:,:,2) = tau .* Jout2(1:L,1:M,1:N);
    
    Jout3 = gather( ifftn( gpuArray(Jout3) ) );    
    JOut(:,:,:,3) = tau .* Jout3(1:L,1:M,1:N);
    
end

% apply op and return local coordinates
JOut = dV * JIn(idx) - JOut(idx);

% clear data
clear JIn; clear tau; clear Jout1; clear Jout2; clear Jout3; clear fJ;




