function [fG] = Generate_Circulant_Nop(G_mn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CIRCULANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_Circulant = tic;
[Gp_mn] = Circulant(G_mn); 
Time_Circulant = toc(tic_Circulant);
fprintf('Time_Circulant   = %dm%ds  \n' ,floor(Time_Circulant/60),int64(mod(Time_Circulant,60)))
%
clear G_mn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FFT - Gmn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_FFT = tic;
[fG] = FFT_Gmn(Gp_mn); 
Time_FFT = toc(tic_FFT);
fprintf('Time_FFT         = %dm%ds  \n' ,floor(Time_FFT/60),int64(mod(Time_FFT,60)))

