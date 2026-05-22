function [Eout,Hout] = field_computation_incident(fG, fK, Jin, dV, f, L ,M ,N, idxI, idxO)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates the incident field Eout and Hout in idxO region
%       due to a current excitation JI in the idxI region
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%


mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*f;


% transform VI in local coordinates to global coordinates
J = zeros(3*L*M*N,1);
J(idxI) = Jin;

if ~isempty(fG)
    [E] = E_field_Nop_comp(J, fG,  dV, omega, eo, zeros(L,M,N,3));
    % get the fields in the desired region idxO
    Eout = E(idxO);
else
    Eout = [];
end

if ~isempty(fK)
    [H] = H_field_Kop_comp(J, fK,  dV, zeros(L,M,N,3));
    % get the fields in the desired region idxO
    Hout = H(idxO);
else
    Hout = [];
end

