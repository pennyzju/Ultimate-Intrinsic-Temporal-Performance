function [Eout,Hout] = evalDyads(Jin,Ocoord,Icoord,ko,dV)
% _________________________________________________________________________
%
%   Dyadic Green Functions
%
% _________________________________________________________________________
%

% % % EM constants
mu = 4*pi*1e-7;
co = 299792458;
% eo = 1/co^2/mu;
%
lambda = 2*pi/ko;
f  = co/lambda;
omega = 2 * pi * f;

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);
Ne = size(Jin,2);

% allocate space for output field
Hx = zeros(No,Ne); % x coomponent, y component, z component
Ex = zeros(No,Ne); % x coomponent, y component, z component
Hy = zeros(No,Ne); % x coomponent, y component, z component
Ey = zeros(No,Ne); % x coomponent, y component, z component
Hz = zeros(No,Ne); % x coomponent, y component, z component
Ez = zeros(No,Ne); % x coomponent, y component, z component

% reshape Jin into 3 components
Jin = reshape(Jin,Ni,3,Ne); % x coomponent, y component, z component

J1 = zeros(1,Ne);
J2 = zeros(1,Ne);
J3 = zeros(1,Ne);

% loop on the Icoord
for jj = 1:Ni
    
    % Get distance vectors
    X = Ocoord(:,1) - Icoord(jj,1);
    Y = Ocoord(:,2) - Icoord(jj,2);
    Z = Ocoord(:,3) - Icoord(jj,3);
         
    % distance 3D
    R2 = X.*X + Y.*Y + Z.*Z;
    R = sqrt(R2);
    R3 = R.*R2;
    
    % For E DGF
    % compute chi
    % chi = 1j*eta*dV*exp(-1j*ko*R)./(4*pi*ko*R3);
    chi = 1j*omega*mu*dV*exp(-1j*ko*R)./(4*pi*ko*ko*R3);
    
    % compute P and Q
    P = 1j*ko*R + 1;
    Q = ko*ko*R2 - P;
    P = (Q - 2*P)./(R2);
    Q = chi.*Q;
    P = chi.*P;
    
    clear chi; clear R2;
    
    % For H DGF
    % precompute value
    const = 1j*ko*R;
    
    % compute mult
    mult = exp(-const)./(4*pi*R3);
    mult = dV*mult.*(const+1);
    
    clear const; clear R; clear R3;
    
    % compute the Solution

    % get components of Jin
    J1(1,:) = Jin(jj,1,:);
    J2(1,:) = Jin(jj,2,:);
    J3(1,:) = Jin(jj,3,:);

    % E field    
    
    % x component
    Ex(:,:) = Ex + (P.*X.*X - Q)*J1 + (P.*(X.*Y))*J2 + (P.*(X.*Z))*J3;
    % y component
    Ey(:,:) = Ey + (P.*(X.*Y))*J1 + (P.*(Y.*Y) - Q)*J2 + (P.*(Y.*Z))*J3;
    % z component
    Ez(:,:) = Ez + (P.*(X.*Z))*J1 + (P.*(Y.*Z))*J2 + (P.*(Z.*Z) - Q)*J3;
    
    
    % H field
    
    X = mult.*X;
    Y = mult.*Y;
    Z = mult.*Z;
    
    % x component
    Hx(:,:) = Hx + Z*J2 - Y*J3;
    % y component
    Hy(:,:) = Hy - Z*J1 + X*J3;
    % z component
    Hz(:,:) = Hz + Y*J1 - X*J2;
    
    
end

Eout = [Ex; Ey; Ez];
clear Ex; clear Ey; clear Ez;
Hout = [Hx; Hy; Hz];
clear Hx; clear Hy; clear Hz;

