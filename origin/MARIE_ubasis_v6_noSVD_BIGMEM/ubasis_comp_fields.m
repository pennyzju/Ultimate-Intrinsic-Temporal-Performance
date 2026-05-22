

function ubasis_comp_fields


clear all;
load ubasis_options;

filename = sprintf('%s/object_def.mat',out_dir);
load(filename);

% load basis solution
Out_file = sprintf('%s/BASIS_Jsol.mat', out_dir);
load(Out_file);
Nexc = size(Jsol,2);

log_filename = sprintf('%s/LOG_comp_fields.txt',out_dir);
fid = fopen(log_filename,'w');


% -------------------------------------------------------------------------
%            DATA
% -------------------------------------------------------------------------

[L,M,N] = size(epsilon_r);
nD = L*M*N;
nS = length(idxS);
nI = length(idxI);
nN = length(idxN);
dV = dx^3;

idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid
idxN3 = [idxN; nD+idxN; 2*nD+idxN]; % the vector of noise related positions in 3D grid


fprintf(fid, '\n\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n');
fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n     # DOFS:                %d', 3*L*M*N);
fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
fprintf(fid, '\n     # DOFS in Coil Shell:  %d', 3*length(idxI));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');
    
type(log_filename);

% -------------------------------------------------------------------------
%                  EM constants
% -------------------------------------------------------------------------

gamma = 42.576e6;

Tesla = f/gamma;
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * f;
lambda  = co/f;
ko = 2*pi/lambda;
omega_mu = omega*mu;
eta =  3.767303134617706e+002; % Free-space impedance

e_r = epsilon_r - 1j*sigma_e/(eo*omega);


% -------------------------------------------------------------------------
%            get the Circulants for the radiated fields
% -------------------------------------------------------------------------

ttot = tic;

% operators with LMN original size

% get the required circulant for the EF operator
[fN] = retrieve_fG(L,M,N,dx,f);

timecirc = toc(ttot);
fprintf(fid,'\n\n Circulants generated, time %f', timecirc);
type(log_filename);

% -------------------------------------------------------------------------
% Apply the direct operator on the FN case
% -------------------------------------------------------------------------

% generate the random excitation matrix, and initialize the fields
fprintf(fid, '\n\n N Operator starting with %d excitations \n', Nexc);
type(log_filename);

% GPU use
dev = gpuDevice(2);
reset(dev);
% send data to GPU mem
% FNgpu = gpuArray(fN);

% load basis
Out_file = sprintf('%s/BASIS_Einc.mat', out_dir);
load(Out_file);

% allocate space
ve = zeros(L,M,N,3);
vj = zeros(L,M,N,3);
Nexc = size(Ebasis,2);

Enoise = zeros(3*nN,Nexc);

t1 = tic;
t2 = tic;
for ii = 1:Nexc
    
    ve(idxS3) = Ebasis(:,ii); % transfer to global vars
    vj(idxS3) = Jsol(:,ii); % transfer to global vars
    
    % compute total field radiated
    [Ein] = E_field_Nop_compGPU_bastien(vj, fN, dV, omega, eo, ve);
    
    % transfer to local vars
    Ebasis(:,ii) = Ein(idxS3);
    Enoise(:,ii) = Ein(idxN3);
    
    if ((ii/100) == floor(ii/100))
        fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
        type(log_filename);
        t2 = tic;
    end
    
end


% FOR HYPERTHERMIA ONLY!!!
save -v7.3 BASIS_Etot.mat Ebasis;

% clear FNgpu;
reset(dev);
clear Ein; clear ve; clear vj;
clear fN;

timeinc = toc(t1);
fprintf(fid,'\n\n Total fields generated generated, time %f', timeinc);
type(log_filename);

% -------------------------------------------------------------------------
%            compute the loss and gsar matrices
% -------------------------------------------------------------------------

t1 = tic;

Ebasis = reshape(Ebasis,nS,3,Nexc); % reshape in x y and z components
Enoise = reshape(Enoise,nN,3,Nexc); % reshape in x y and z components

Sig = sigma_e(idxS);
Rho = rho(idxS);

SigN = sigma_e(idxN);
RhoN = rho(idxN);

% tic;
% SIGDIAG = spdiags(Sig,0,nS,nS);
% LOSS = Ex.'*SIGDIAG*conj(Ex) + Ey.'*SIGDIAG*conj(Ey) + Ez.'*SIGDIAG*conj(Ez);
% toc
% 
% tic
% RHODIAG = spdiags(Sig./Rho,0,nS,nS);
% GSAR = Ex.'*RHODIAG*conj(Ex) + Ey.'*RHODIAG*conj(Ey) + Ez.'*RHODIAG*conj(Ez);
% toc


SIGDIAG = spdiags(Sig,0,nS,nS) * pvol;
LOSS = (squeeze(Ebasis(:,1,:)))'*(SIGDIAG*(squeeze(Ebasis(:,1,:)))) + (squeeze(Ebasis(:,2,:)))'*(SIGDIAG*(squeeze(Ebasis(:,2,:)))) + (squeeze(Ebasis(:,3,:)))'*(SIGDIAG*(squeeze(Ebasis(:,3,:))));
LOSS = LOSS.';  % LOSS is in W

RHODIAG = spdiags(Sig./Rho,0,nS,nS);
GSAR = (squeeze(Ebasis(:,1,:)))'*(RHODIAG*(squeeze(Ebasis(:,1,:)))) + (squeeze(Ebasis(:,2,:)))'*(RHODIAG*(squeeze(Ebasis(:,2,:)))) + (squeeze(Ebasis(:,3,:)))'*(RHODIAG*(squeeze(Ebasis(:,3,:))));

npix = size(idxS,1);
GSAR = GSAR.' / (2*npix);  % SAR is in W/kg averaged over the whole volume (hence the factor npix)


clear Ebasis;

t2 = toc(t1);
fprintf(fid, '\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
type(log_filename);

Xout_file = sprintf('%s/BASIS_LOSS.mat', out_dir);
save(Xout_file, 'LOSS', '-v7.3');
Xout_file = sprintf('%s/BASIS_GSAR.mat', out_dir);
save(Xout_file, 'GSAR', '-v7.3');

t1 = tic;

SIGDIAG = spdiags(SigN,0,nN,nN) * pvol;
LOSS = (squeeze(Enoise(:,1,:)))'*(SIGDIAG*(squeeze(Enoise(:,1,:)))) + (squeeze(Enoise(:,2,:)))'*(SIGDIAG*(squeeze(Enoise(:,2,:)))) + (squeeze(Enoise(:,3,:)))'*(SIGDIAG*(squeeze(Enoise(:,3,:))));
LOSS = LOSS.';  % LOSS is in W

RHODIAG = spdiags(SigN./RhoN,0,nN,nN);
GSAR = (squeeze(Enoise(:,1,:)))'*(RHODIAG*(squeeze(Enoise(:,1,:)))) + (squeeze(Enoise(:,2,:)))'*(RHODIAG*(squeeze(Enoise(:,2,:)))) + (squeeze(Enoise(:,3,:)))'*(RHODIAG*(squeeze(Enoise(:,3,:))));

npix = size(idxN,1);
GSAR = GSAR.' / (2*npix);  % SAR is in W/kg averaged over the whole volume (hence the factor npix)

t2 = toc(t1);
fprintf(fid, '\n\n LOSS and GSAR matrices computed. Elapsed time %f \n', Nexc, t2);
type(log_filename);

Xout_file = sprintf('%s/BASIS_LOSS_N.mat', out_dir);
save(Xout_file, 'LOSS', '-v7.3');
Xout_file = sprintf('%s/BASIS_GSAR_N.mat', out_dir);
save(Xout_file, 'GSAR', '-v7.3');

% -------------------------------------------------------------------------
%            get the Circulants for the radiated fields
% -------------------------------------------------------------------------

ttot = tic;

% operators with LMN original size

% get the required circulant for the MF operator
[fK] = retrieve_fK(L,M,N,dx,f);

timecirc = toc(ttot);
fprintf(fid,'\n\n Circulants generated, time %f', timecirc);
type(log_filename);

% -------------------------------------------------------------------------
% Apply the direct operator on the FK case
% -------------------------------------------------------------------------


% generate the random excitation matrix, and initialize the fields
fprintf(fid, '\n\n K Operator starting with %d excitations \n', Nexc);
type(log_filename);

% GPU use
dev = gpuDevice(1);
reset(dev);
% send data to GPU mem
% FKgpu = gpuArray(fK);

% load basis
Out_file = sprintf('%s/BASIS_Hinc.mat', out_dir);
load(Out_file);

% allocate space
vm = zeros(L,M,N,3);
vj = zeros(L,M,N,3);
Nexc = size(Hbasis,2);

t1 = tic;
t2 = tic;

for ii = 1:Nexc
    
    vm(idxS3) = Hbasis(:,ii); % transfer to global vars
    vj(idxS3) = Jsol(:,ii); % transfer to global vars
    
    % compute total field radiated
    % incident field is zero (0*ve)
    [Hin] = H_field_Kop_compGPU_bastien(vj, fK, dV, vm);
    
    % transfer to local vars
    Hbasis(:,ii) = mu*Hin(idxS3);
    
    if ((ii/500) == floor(ii/500))
        fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
        type(log_filename);
        t2 = tic;
    end
    
end

% clear FKgpu; 
reset(dev);
clear vj; clear ve; clear Hin;
clear fK;

timeinc = toc(t1);
fprintf(fid,'\n\n Total fields generated generated, time %f', timeinc);
type(log_filename);

% -------------------------------------------------------------------------
%            compute the B1p and B1m matrices
% -------------------------------------------------------------------------

t1 = tic;

Hbasis = reshape(Hbasis,nS,3,Nexc); % reshape in x y and z components

B1p = ( squeeze(Hbasis(:,1,:)) + 1j*squeeze(Hbasis(:,2,:)) ) / sqrt(2);
Xout_file = sprintf('%s/BASIS_B1p.mat', out_dir);
save(Xout_file, 'B1p', '-v7.3');

B1m = ( squeeze(Hbasis(:,1,:)) - 1j*squeeze(Hbasis(:,2,:)) ) / sqrt(2);
Xout_file = sprintf('%s/BASIS_B1m.mat', out_dir);
save(Xout_file, 'B1m', '-v7.3');

clear Hbasis;
clear B1p; clear B1m;

t2 = toc(t1);
fprintf(fid, '\n\n B1p and B1m matrices computed. Elapsed time %f \n', t2);
type(log_filename);

fclose(fid);









