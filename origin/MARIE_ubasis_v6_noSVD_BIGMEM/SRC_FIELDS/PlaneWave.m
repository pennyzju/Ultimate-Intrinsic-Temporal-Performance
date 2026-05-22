function [Eexc,Hexc] = PlaneWave(dx,dy,dz,r,ko,polarization,direction)

%
[L, M, N, ~] = size(r);

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
eta = sqrt(mu/eo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Excitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define excitation - a plane wave
if strcmp(direction,'x') 
    Eexc = dy*dz* exp(-1i*ko* r(:,:,:,1) ) *  (exp(-1i*ko*dx/2) - exp(1i*ko*dx/2)) / (-1i*ko)  ;
elseif strcmp(direction,'y') 
    Eexc = dx*dz* exp(-1i*ko* r(:,:,:,2) ) *  (exp(-1i*ko*dy/2) - exp(1i*ko*dy/2)) / (-1i*ko)  ;
elseif strcmp(direction,'z') 
    Eexc = dx*dy* exp(-1i*ko* r(:,:,:,3) ) *  (exp(-1i*ko*dz/2) - exp(1i*ko*dz/2)) / (-1i*ko)  ;
end

if strcmp(polarization,'x') 
    Hexc = [zeros(L*M*N,1) ; Eexc(:) ; zeros(L*M*N,1) ]./eta;
    Eexc = [Eexc(:) ; zeros(2*L*M*N,1) ];    
elseif strcmp(polarization,'y')
    Eexc = [zeros(L*M*N,1) ; Eexc(:) ; zeros(L*M*N,1) ];
elseif strcmp(polarization,'z')
    Eexc = [zeros(2*L*M*N,1) ; Eexc(:) ];
end

%
Eexc = reshape(Eexc , L, M, N, 3);
Hexc = reshape(Hexc , L, M, N, 3);