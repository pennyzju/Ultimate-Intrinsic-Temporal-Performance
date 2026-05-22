function [epsilon_r,sigma_e,rho] = homogensphere(r,Rc,Radius,epsilon,sigma,density)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates a homogeneous sphere
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%

[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
sigma_e = zeros(L,M,N);
rho = zeros(L,M,N);

% define sphere
sphere = @(r)( (r(:,:,:,1) - Rc(1,1) ).^2 + ( r(:,:,:,2) - Rc(2,1) ).^2 + ( r(:,:,:,3) - Rc(3,1) ).^2 < Radius^2) ;

pointsphere = sphere(r); % ones for domain elements in the shell
idx = find(pointsphere(:)); % get indexes of elements

epsilon_r(idx) = epsilon; 
sigma_e(idx) = sigma; 
rho(idx) = density; 
