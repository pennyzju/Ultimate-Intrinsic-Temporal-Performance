function [JOut] = mv_op_fft_Nop_1(JIn0, fG, e_r, dV, transp_flag, idx)
%
% Function that applies the JVIE I core operation
%   Computational Prototyping Group, RLE at MIT
%


% fft dimensions
LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);

[L, M, N] = size(e_r);

% allocate space
JIn = zeros(L, M, N, 3);
JOut = zeros(L, M, N, 3);

% translate from local to global coordinates
JIn(idx) = JIn0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(transp_flag,'transp')       % y = A'*x
    
    chi = conj((e_r - 1.0));

       
    % apply fft and mv-op
    fJ = fftn(conj( chi .* JIn(:,:,:,1)) , [LfG, MfG, NfG]);
    Jout1 = fG(:,:,:,1) .* fJ;
    Jout2 = fG(:,:,:,2) .* fJ;
    Jout3 = fG(:,:,:,3) .* fJ;
    
    fJ = fftn(conj( chi .* JIn(:,:,:,2)) , [LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,5) .* fJ;
    
    fJ = fftn(conj( chi .* JIn(:,:,:,3)) , [LfG, MfG, NfG]);
    Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fG(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = dV * conj(e_r) .* JIn(:,:,:,1) - conj( Jout1(1:L,1:M,1:N));
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = dV * conj(e_r) .* JIn(:,:,:,2) - conj( Jout2(1:L,1:M,1:N));
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = dV * conj(e_r) .* JIn(:,:,:,3) - conj( Jout3(1:L,1:M,1:N));
    


elseif strcmp(transp_flag,'notransp') % y = A*x
    
    chi = e_r - 1.0 ;
        
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
    JOut(:,:,:,1) = dV * e_r .* JIn(:,:,:,1) - chi .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = dV * e_r .* JIn(:,:,:,2) - chi .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = dV * e_r .* JIn(:,:,:,3) - chi .* Jout3(1:L,1:M,1:N);
    

end

% return local coordinates
JOut = JOut(idx);