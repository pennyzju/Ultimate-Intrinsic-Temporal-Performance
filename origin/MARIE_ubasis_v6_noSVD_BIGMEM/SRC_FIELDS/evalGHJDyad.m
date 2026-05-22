function [Hout] = evalGHJDyad(Jin,Ocoord,Icoord,ko,dV)
% _________________________________________________________________________
%
%   Dyadic Green Function
%
% _________________________________________________________________________
%

% % % EM constants
% % mu = 4*pi*1e-7;
% % co = 299792458;
% % eo = 1/co^2/mu;
% % %
% % lambda = 2*pi/ko;
% % f  = co/lambda;
% % omega = 2 * pi * f;


% % Free-space impedance
% eta =  3.767303134617706e+002; 

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);

% allocate space for output field
Hout = zeros(No,3); % x coomponent, y component, z component

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
clear R2;

% precompute value
const = 1j*ko*R;

% compute mult
mult = exp(-const)./(4*pi*R3);
mult = dV*mult.*(const+1);


clear const; clear R; clear R2; clear R3;

% compute the Solution

X = mult.*X;
Y = mult.*Y;
Z = mult.*Z;

% x component
Hout(:,1) = Z*Jin(:,2) - Y*Jin(:,3);

% y component
Hout(:,2) = -Z*Jin(:,1) + X*Jin(:,3);

% z component
Hout(:,3) = Y*Jin(:,1) - X*Jin(:,2);


