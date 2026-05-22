function [Eout] = evalGEJDyad(Jin,Ocoord,Icoord,ko,dV)
% _________________________________________________________________________
%
%   Dyadic Green Function
%
% _________________________________________________________________________
%

% EM constants
mu = 4*pi*1e-7;
co = 299792458;
% eo = 1/co^2/mu;
%
lambda = 2*pi/ko;
f  = co/lambda;
omega = 2 * pi * f;

% Free-space impedance
% eta =  3.767303134617706e+002; 

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);

% allocate space for output field
Eout = zeros(No,3); % x coomponent, y component, z component

% reshape Jin into 3 components
Jin = reshape(Jin,Ni,3); % x coomponent, y component, z component

% Get distance vectors
X = zeros(No,Ni); Y = zeros(No,Ni); Z = zeros(No,Ni);
for jj = 1:Ni
    X(:,jj) = Ocoord(:,1) - Icoord(jj,1);
    Y(:,jj) = Ocoord(:,2) - Icoord(jj,2);
    Z(:,jj) = Ocoord(:,3) - Icoord(jj,3);
end

% distance 3D
R2 = X.*X + Y.*Y + Z.*Z;
R = sqrt(R2);
R3 = R.*R2;

% compute chi
chi = 1j*omega*mu*dV*exp(-1j*ko*R)./(4*pi*ko*ko*R3);
% chi = 1j*eta*dV*exp(-1j*ko*R)./(4*pi*ko*R3);

% compute P and Q
P = 1j*ko*R + 1;
Q = ko*ko*R2 - P;
P = (Q - 2*P)./(R2);
Q = chi.*Q;
P = chi.*P;

clear chi; clear R; clear R2; clear R3;

% compute the Solution
 
% x component
Eout(:,1) = (P.*X.*X - Q)*Jin(:,1) + (P.*(X.*Y))*Jin(:,2) + (P.*(X.*Z))*Jin(:,3);

% y component
Eout(:,2) = (P.*(X.*Y))*Jin(:,1) + (P.*(Y.*Y) - Q)*Jin(:,2) + (P.*(Y.*Z))*Jin(:,3);

% z component
Eout(:,3) = (P.*(X.*Z))*Jin(:,1) + (P.*(Y.*Z))*Jin(:,2) + (P.*(Z.*Z) - Q)*Jin(:,3);

