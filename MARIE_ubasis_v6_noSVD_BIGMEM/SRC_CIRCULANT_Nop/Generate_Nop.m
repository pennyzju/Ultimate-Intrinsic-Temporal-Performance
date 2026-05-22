function [G_mn] = Generate_Nop(r,ko,dx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Define grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, M, N, ~] = size(r);
dy =dx; dz = dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic_Assembly = tic;

% Set order of quadrature for each region
Np_DEMCEM    = 20;  
Np_1D_near   = 20;
Np_1D_medium = 10;
Np_1D_far    = 3;
Np_GL = [Np_DEMCEM ; Np_1D_near ; Np_1D_medium ; Np_1D_far];

% Set region of medium distances cells for higher order quadrature
n_medium = 5; 

% Get the 6 faces of  the cube
[R_faces] = Cube_faces(dx,dy,dz);

% Assembly
[G_mn] = Assembly_Nop(L,M,N,Np_GL,n_medium,ko,dx,dy,dz,R_faces);
%

Time_Assembly = toc(tic_Assembly);
fprintf('------------------------  \n')
fprintf('Time_Assembly    = %dm%ds  \n' ,floor(Time_Assembly/60),int64(mod(Time_Assembly,60)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Close parallelization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%