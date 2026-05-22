function [Uin,Sin,Xin,Pin,xds,yds,zds] = MRGF_incident_ROM(e_r, r, f, Iregion, tol, outfile, nexc)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates the MRGF for a given domain and eps
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

fid = 1;
% the output file will be in the ..\DATA_MRGF folder
outpath = strcat('.\DATA_MRGF\', outfile);


% -------------------------------------------------------------------------
%                  Prepare the cases to Solve
% -------------------------------------------------------------------------

[L,M,N] = size(e_r); % get dimensions of the problem


idxS = find(abs(e_r(:))-1);

if isempty(Iregion)
    idxI = find((abs(e_r(:))-1 )== 0); % get indexes of air elements
else
    idxI = Iregion(r); % get indexes of elements in vector of I domain
end

nD = L*M*N;
nS = length(idxS);
nI = length(idxI);

idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid

dx = r(2,1,1,1) - r(1,1,1,1);
dV = dx.^3;

% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * f;
lambda  = co/f;
ko = 2*pi/lambda;
omega_mu = omega*mu;
eta =  3.767303134617706e+002; % Free-space impedance


% % -------------------------------------------------------------------------
% %                  Generate the Circulant
% % -------------------------------------------------------------------------
% 
% % Multiples of 16
% LfG = 4*(ceil(L/4));
% MfG = 4*(ceil(M/4));
% NfG = 4*(ceil(N/4));
% 
% [fG] = retrieve_fG(LfG,MfG,NfG,dx,f);


% -------------------------------------------------------------------------
%                  Save data
% -------------------------------------------------------------------------

DATA_file = sprintf('%s_MRGF_DATA.mat', outpath);
save(DATA_file, 'fG', 'r', 'e_r', 'idxI', 'f', '-v7.3');


% -------------------------------------------------------------------------
% Plot geometry
% -------------------------------------------------------------------------

% xd = r(:,:,:,1);
% yd = r(:,:,:,2);
% zd = r(:,:,:,3);
% 
% plot3(xd(idxS),yd(idxS),zd(idxS),'b.');
% hold on;
% plot3(xd(idxI),yd(idxI),zd(idxI),'g.');
% pause;


% -------------------------------------------------------------------------
% Generate Incident Field
% -------------------------------------------------------------------------

tic

% AMIT's DDA N circulant
[fG] = dda_Ncirculant(r, dx, ko);

fprintf(fid,'\n\n Circulant generated, computing the incident fields for %d excitations', nexc);

Ee = zeros(3*nS,nexc);
ve = zeros(L,M,N,3);

t1 = clock;

for ii = 1:nexc
    
    %     tic;
    Je = rand(3*nI,1) + 1j*rand(3*nI,1);
    ve(idxI3) = Je/norm(Je); % transfer to global vars
    
    [Ein] = dda_Eradiate(ve, fG, ko, L, M ,N);
    Ee(:,ii) = Ein(idxS3);
        
    %         timedyad = toc;
    %     fprintf(fid,'\n iteration %d, time %f', ii, timedyad);
    if ((ii/500) == floor(ii/500))
        t2 = clock;
        fprintf(fid, '\n %d Fields. Elapsed time %f', ii, etime(t2,t1));
    end
    
    
    
end

clear Je; clear ve; clear Ein;

Basis_file = sprintf('%s_MRGF_Ee.mat', outpath);
save(Basis_file, 'Ee', '-v7.3');

timeinc = toc;
fprintf(fid,'\n\n Incident fields generated generated, time %f', timeinc);

% -------------------------------------------------------------------------
% Compute the SVD of the basis
% -------------------------------------------------------------------------

tic;

[Qe, ~, ~] = svd(Ee,'econ');
clear Ee;

Basis_file = sprintf('%s_MRGF_Qe.mat', outpath);
save(Basis_file, 'Qe', '-v7.3');

timesvd = toc;
fprintf(fid,'\n SVD done, time %f\n', timesvd);


% -------------------------------------------------------------------------
%            apply the rest of the RSVD
% -------------------------------------------------------------------------

tic;


% apply again the FIELD_COMPUATION function, in this case the adjoint op
Uin = zeros(size(Qe,2),3*nI); % B is going to be Ui
ve = zeros(L,M,N,3);
for localii = 1:size(Qe,2)
    
    ve(idxS3) = conj(Qe(:,localii)); % transfer to global vars
    
    [Ein] = dda_Eradiate(ve, fG, ko, L, M ,N);
    Uin(localii,:) = (Ein(idxI3)).';

end
clear Qe;
timeinc = toc;
fprintf(fid,'\n Adjoint operator done, time %f', timeinc);

tic
% apply the SVD to generate the basis and singular values
[Uin, Sin, ~] = svd(Uin,'econ');
Sin = sparse(Sin); 

timeinc = toc;
fprintf(fid,'\n SVD done, time %f', timeinc);

Basis_file = sprintf('%s_MRGF_SV.mat', outpath);
save(Basis_file, 'Sin', '-v7.3');

tic

Basis_file = sprintf('%s_MRGF_Qe.mat', outpath);
load(Basis_file);

Uin = Qe*Uin; 
clear Qe;

Basis_file = sprintf('%s_MRGF_Uin.mat', outpath);
save(Basis_file, 'Uin', '-v7.3');


timeinc = toc;
fprintf(fid,'\n Right MVP and store done, time %f', timeinc);

fprintf(fid, '\nRandom SVD operation done with %d excitations\n', size(Uin,2));



% -------------------------------------------------------------------------
% Truncate the singular values
% -------------------------------------------------------------------------

numsv = length(Sin);
for ii = 2:length(Sin)
    if (abs(Sin(ii,ii)) < tol*abs(Sin(1,1))) || (abs(Sin(ii,ii)) < 5e-4)
        numsv = ii;
        break;
    end
end

fprintf(fid, '\n Truncation to %d singular vectors', numsv);


% -------------------------------------------------------------------------
% Apply the DEIM procedure on the incident basis
% -------------------------------------------------------------------------

ta = clock;
% We will use up to 3 times the number of singular values
coeffdeim = min([size(Uin,2), 3*numsv]);
Uin = Uin(:,1:coeffdeim);

[~, Pin] = deim(Uin);

% transform the DEIM points into DEIM voxels: i.e. add the 3 components of
% the selected points

% get S domain coordinates
xs = r(:,:,:,1); xs = xs(idxS);
ys = r(:,:,:,2); ys = ys(idxS);
zs = r(:,:,:,3); zs = zs(idxS);

% find all the voxels with at least one component selected by deim
Pin3D = [Pin(1:nS,:), Pin(nS+1:2*nS,:), Pin(2*nS+1:3*nS,:)];
vec3D = sum(Pin3D,2);
idxD = find(vec3D);
idxD = sort(idxD);

% get deim points coordinates
xds = xs(idxD);
yds = ys(idxD);
zds = zs(idxD);

% 3D idx and build the Pin for the DEIM approach
idxD = [idxD; nS+idxD; 2*nS+idxD];
e = speye(3*nS);
Pin = sparse(e(:, idxD));
ndeim = size(Pin,2);

% Obtain the DEIM extended inverse matrix
Xin= (Pin.'*Uin)\speye(ndeim,ndeim);

% do the actual truncation of the basis
Uin = Uin(:,1:numsv);
Sin = Sin(1:numsv,1:numsv);

% truncate also the DEIM Xdx extended matrix
Xin = Xin(1:numsv,:);

tb = clock;
fprintf(fid, '\n\n DEIM algorithm applied on %d vectors. Elapsed time %f \n', ndeim, etime(tb,ta));

DEIM_file = sprintf('%s_MRGF_DEIM.mat', outpath);
save(DEIM_file, 'Xin', 'Pin', 'xds', 'yds', 'zds', '-v7.3');

IBasis_file = sprintf('%s_MRGF_Utrunc.mat', outpath);
save(IBasis_file, 'Uin', 'Sin', '-v7.3');


% matlabpool close;
