function [fG] = dda_Ncirculant(r, dx, k0)
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Amit's code for generating the circulant
%
% _________________________________________________________________________
% _________________________________________________________________________


% Some EM data
eta =  3.767303134617706e+002; % Free-space impedance

% get grid size
[L, M, N, ~] = size(r);

% compute second grid
x = ((1:2*L-1)-L)*dx; 
y = ((1:2*M-1)-M)*dx;
z = ((1:2*N-1)-N)*dx;
r2 = grid3d(x, y, z);
clear x; clear y; clear z;

a = dx*(3/4/pi)^(1/3); % Radius of equivalent sphere
 
R = dda_abs3d(r2);
unitR = dda_scalar_mult(1./R, r2);
alpha = k0*R;
alpha2 = alpha.^2;
alpha3 = alpha.^3;

fG = zeros(2*L, 2*M, 2*N, 3, 3);

% Calc dyadic Green's function
for q = 1:3
    for p = 1:3
        if p == q
            deltaFact = (alpha2-1-1j*alpha);
        else
            deltaFact = 0;
        end
        G = -1j*k0^2*eta*exp(-1j*alpha)./(4*pi*alpha3).*...
            (deltaFact + unitR(:,:,:,p).*unitR(:,:,:,q).*(3-alpha2+3j*alpha));
        if p==q
            G(L,M,N) = -2j*eta/k0/3*(exp(-1j*k0*a)*(1+1j*k0*a)-1);
        else
            G(L,M,N) = 0;
        end
        fG(:,:,:,p,q) = fftn(G, [2*L, 2*M, 2*N]);
    end
end

