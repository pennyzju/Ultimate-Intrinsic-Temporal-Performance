
function ubasis_comp_inc_basis_bigmask(public_out_dir,out_dir,fileName)

    % load data from file
load(fullfile(public_out_dir, 'ubasis_options.mat'));
load(fullfile(out_dir,fileName));

log_filename = sprintf('%s/LOG_comp_inc_fiels.txt',public_out_dir);
fid = fopen(log_filename,'w');
tini = tic;

gamma = 42.576e6;

% -------------------------------------------------------------------------
%            DATA
% -------------------------------------------------------------------------

filename = sprintf('%s/object_def_bigmask.mat',out_dir);
load(filename);
[L,M,N] = size(rotated_epsilon_r);

nD = L*M*N;
nS = length(idxS); %需要保存的数据点
nI = length(idxI);
dV = dx^3;

idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid 偶极子云的位置
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid 人头的位置

fprintf(fid, '\n\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n');
fprintf(fid, '\n     Domain:                %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n     Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n     # DOFS in Domain:      %d', 3*L*M*N);
fprintf(fid, '\n     # DOFS in Scatterer:   %d', 3*length(idxS));
fprintf(fid, '\n     # DOFS in Dipoles:     %d', 3*length(idxI));
fprintf(fid, '\n     Frequency:             %.2fMHz', f/1e6);
fprintf(fid, '\n     # Random Excitations:  %d', Nexc);
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');
    
type(log_filename);


% -------------------------------------------------------------------------
%                  EM constants
% -------------------------------------------------------------------------


Tesla = f/gamma;
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;

omega = 2 * pi * f;
lambda  = co/f;
ko = 2*pi/lambda;
omega_mu = omega*mu;
eta =  3.767303134617706e+002; % Free-space impedance

ce = 1j*omega*eo;
cm = 1j*omega*mu;

% e_r = epsilon_r - 1j*sigma_e/(eo*omega);%计算复数介电常数


ttot = tic;

% get the required circulant for the EF operator
[fN] = retrieve_fG(L,M,N,dx,f);

timecirc = toc;
fprintf(fid,'\n\n Circulant fN generated, time %f', timecirc);%fprintf 函数记录并写入一条日志信息到指定文件中。

fprintf('\n\n Circulant fN generated, time %f', timecirc);
% type(log_filename);


% -------------------------------------------------------------------------
%            Generate the random excitation bewteen -1 and 1
% -------------------------------------------------------------------------
% 设置随机数生成器的种子
seed = 42;  % 你可以选择任意整数作为种子
rng(seed);


% 检查 Je 是否存在
Out_file = sprintf('%s/Je.mat', public_out_dir);

if exist(Out_file, 'file') == 2
    % 如果 Je.mat 存在，加载它
    load(Out_file, 'Je');
else
    % 如果 Je.mat 不存在，创建 Je 并保存
    Je = 2*(rand(3*nI,Nexc)-0.5) + 1j*2*(rand(3*nI,Nexc)-0.5);
    save(Out_file , 'Je', '-v7.3');
end


% -------------------------------------------------------------------------
% Apply the direct operator on the FN case
% -------------------------------------------------------------------------


% generate the random excitation matrix, and initialize the fields
fprintf(fid, '\n\n N Operator starting with %d excitations \n', Nexc);
fprintf('\n\n N Operator starting with %d excitations \n', Nexc);
% type(log_filename);
fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
% GPU use
dev = gpuDevice(1);
reset(dev);
% send data to GPU mem
% FNgpu = gpuArray(fN);

% allocate space
ve = zeros(L,M,N,3);

Ebasis = ones(3*nS,2*Nexc) + 1j*ones(3*nS,2*Nexc);  % store both the E and H dipoles (2*Nexc)
Hbasis = ones(3*nS,2*Nexc) + 1j*ones(3*nS,2*Nexc);
fprintf('当前系统时间: %s\n', datestr(datetime('now'), 'ddd mmm dd HH:MM:SS yyyy'));
t1 = tic;
t2 = tic;
for ii = 1:Nexc
    
     ve(idxI3) = Je(:,ii)/norm(Je(:,ii)); % transfer to global vars and normalize
    
    % compute total field radiated
    % incident field is zero (0*ve)
    NN_je = E_field_Nop_compGPU_bastien(ve, fN, dV, omega, eo, 0*ve);  % 1/ce*N
    
    % transfer to local vars
    Ebasis( :,2*(ii-1)+1 ) =  NN_je(idxS3);  % E dipole contribution to the E field    
    Hbasis( :,2*(ii-1)+2 ) =  ce/cm*NN_je(idxS3);  % M dipole contribution to the H field
    
    if 1 % ((ii/100) == floor(ii/100))
        fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
        fprintf('\n %d Fields. Elapsed time %g', ii, toc(t2));
        %type(log_filename);
        t2 = tic;
    end
    
end

% clear FNgpu; reset(dev);
clear NN_je;
clear fN;

timeinc = toc(t1);
fprintf(fid,'\n\n Incident fields generated generated, time %f', timeinc);
fprintf('\n\n Incident fields generated generated, time %f', timeinc);
% type(log_filename);


% -------------------------------------------------------------------------
% Apply the direct operator on the FK case
% -------------------------------------------------------------------------

ttot = tic;

% get the required circulant for the MF operator
[fK] = retrieve_fK(L,M,N,dx,f);

timecirc = toc;
fprintf(fid,'\n\n Circulant fK generated, time %f', timecirc);
fprintf('\n\n Circulant fK generated, time %f', timecirc);
% type(log_filename);


% generate the random excitation matrix, and initialize the fields
fprintf(fid, '\n\n K Operator starting with %d excitations \n', Nexc);
fprintf('\n\n K Operator starting with %d excitations \n', Nexc);
% type(log_filename);

% GPU use
dev = gpuDevice(1);
reset(dev);
% send data to GPU mem
% FKgpu = gpuArray(fK);

t1 = tic;
t2 = tic;

for ii = 1:Nexc
    
    ve(idxI3) = Je(:,ii)/norm(Je(:,ii)); % transfer to global vars and normalize
    
    % compute total field radiated
    % incident field is zero (0*ve)
    KK_je = H_field_Kop_compGPU_bastien(ve, fK, dV, 0*ve);  % K
    
    % transfer to local vars
    Hbasis( :,2*(ii-1)+1 ) =  KK_je(idxS3);  % E dipole contribution to the H field
    Ebasis( :,2*(ii-1)+2 ) = -KK_je(idxS3);  % M dipole contribution to the E field   
    
    if 1  % ((ii/100) == floor(ii/100))
        fprintf(fid, '\n %d Fields. Elapsed time %g', ii, toc(t2));
        fprintf('\n %d Fields. Elapsed time %g', ii, toc(t2));
        % type(log_filename);
        t2 = tic;
    end
    
end

% clear FKgpu; reset(dev);
clear Je; clear ve; clear KK_je;
clear fK;

timeinc = toc(t1);
fprintf(fid,'\n\n Incident fields generated generated, time %f', timeinc);
fprintf('\n\n Incident fields generated generated, time %f', timeinc);
% type(log_filename);


% Check if the 'bigmask' directory exists
if ~exist(fullfile(out_dir, 'bigmask'), 'dir')
    % If it doesn't exist, create it
    mkdir(fullfile(out_dir, 'bigmask'));
end
Out_file = sprintf('%s/bigmask/BASIS_Einc.mat', out_dir);
save(Out_file , 'Ebasis', '-v7.3');

Out_file = sprintf('%s/bigmask/BASIS_Hinc.mat', out_dir);
save(Out_file , 'Hbasis', '-v7.3');




