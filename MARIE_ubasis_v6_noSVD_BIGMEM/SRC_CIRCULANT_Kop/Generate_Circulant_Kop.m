function [fK] = Generate_Circulant_Kop(Kop_mn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CIRCULANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_Circulant = tic;
[CK_mn] = Circulant_Kop(Kop_mn); 
Time_Circulant = toc(tic_Circulant);
fprintf('Time_Circulant   = %dm%ds  \n' ,floor(Time_Circulant/60),int64(mod(Time_Circulant,60)))
%
clear Kop_mn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FFT - Gmn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_FFT = tic;
[fK] = FFT_CK(CK_mn); 
Time_FFT = toc(tic_FFT);
fprintf('Time_FFT         = %dm%ds  \n' ,floor(Time_FFT/60),int64(mod(Time_FFT,60)))



