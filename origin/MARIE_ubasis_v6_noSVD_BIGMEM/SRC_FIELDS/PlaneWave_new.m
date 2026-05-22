function [Eexc,Hexc] = PlaneWave_new(r,k,omega_mu,polarization)
%%
dx = r(2,1,1,1) - r(1,1,1,1);
dy = r(1,2,1,2) - r(1,1,1,2);
dz = r(1,1,2,3) - r(1,1,1,3);
% k = (kx,ky,kz)
kx = k(1);
ky = k(2);
kz = k(3);
%
[L, M, N, ~] = size(r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Excitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if kx ~= 0
    Ex =  exp(-1i*kx* r(:,:,:,1) ) *  (exp(-1i*kx*dx/2) - exp(1i*kx*dx/2)) / (-1i*kx)  ;
else
    Ex = dx * ones(L,M,N);
end

if ky ~= 0
    Ey = exp(-1i*ky* r(:,:,:,2) ) *  (exp(-1i*ky*dy/2) - exp(1i*ky*dy/2)) / (-1i*ky)  ;
else
    Ey = dy * ones(L,M,N);
end

if kz ~= 0
    Ez = exp(-1i*kz* r(:,:,:,3) ) *  (exp(-1i*kz*dz/2) - exp(1i*kz*dz/2)) / (-1i*kz)  ;
else
    Ez = dz * ones(L,M,N);
end
%
E = Ex .* Ey .* Ez;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(polarization,'x')

    Eexc = [E(:) ; zeros(2*L*M*N,1) ];
    %
    Hexc =  [zeros(L*M*N,1) ; (kz /omega_mu) * E(:) ; -(ky /omega_mu) * E(:) ];

elseif strcmp(polarization,'y')

    Eexc = [zeros(L*M*N,1) ; E(:) ; zeros(L*M*N,1) ];
    %
    Hexc =  [-(kz /omega_mu) * E(:) ; zeros(L*M*N,1) ;  (kx /omega_mu) * E(:) ];

elseif strcmp(polarization,'z')

    Eexc = [zeros(2*L*M*N,1) ; E(:) ];
    %
    Hexc =  [(ky /omega_mu) * E(:) ; -(kx /omega_mu) * E(:) : zeros(L*M*N,1) ];

end
Eexc = reshape(Eexc , L, M, N, 3);
Hexc = reshape(Hexc , L, M, N, 3);