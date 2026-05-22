function [I_ST] = DEMCEM_ST(Np_1D,ko,dx,dy)
%%

% GEOMETRY
r1 = [0, 0,0]';
r2 = [dx,0,0]';
r3 = [dx,dy,0]';
r4 = [0 ,dy,0]';

[w,z] = Gauss_1D_Nop(Np_1D);
%
N_theta = Np_1D;
N_psi = Np_1D;
%
[ w_theta , z_theta ] = Gauss_1D_Nop ( N_theta );
[ w_psi , z_psi ] = Gauss_1D_Nop ( N_psi );

% Evaluate ST

[I_T1_T1] = mexDIRECT_WS_ST_const(r1,r2,r3,ko,Np_1D,z,w);
[I_T1_T2] = mexDIRECT_WS_EA_const(r1,r3,r2,r4, ko, N_theta, N_psi , w_theta, z_theta, w_psi, z_psi);
[I_T2_T1] = mexDIRECT_WS_EA_const(r1,r3,r4,r2, ko, N_theta, N_psi , w_theta, z_theta, w_psi, z_psi);
[I_T2_T2] = mexDIRECT_WS_ST_const(r1,r3,r4,ko,Np_1D,z,w);

% [I_T1_T1] = DIRECT_WS_ST_RWG(r1,r2,r3,Np_1D);
% [I_T1_T2] = DIRECT_WS_EA_RWG(r1,r3,r2,r4,Np_1D);
% [I_T2_T1] = DIRECT_WS_EA_RWG(r1,r3,r4,r2,Np_1D);
% [I_T2_T2] = DIRECT_WS_ST_RWG(r1,r3,r4,Np_1D);

I_ST = I_T1_T1 + I_T1_T2 + I_T2_T1 + I_T2_T2;
