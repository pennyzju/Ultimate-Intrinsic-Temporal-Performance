function [G_mn] = Assembly_new(L,M,N,Np_GL,n_medium)
%%
addpath('../SRC_CIRCULANT_N/GL_CUBATURE')
addpath('../SRC_CIRCULANT_N/DEMCEM')
addpath('../SRC_CIRCULANT_N/DEMCEM/mexDIRECT_WS_const')

% Define global parameters
global ko
global dx dy dz
global R_faces
% Set order of quadrature for each region
% Np_GL = [Np_DEMCEM ; Np_1D_near ; Np_1D_medium ; Np_1D_far];
Np_DEMCEM    = Np_GL(1);
Np_1D_near   = Np_GL(2);
Np_1D_medium = Np_GL(3);
Np_1D_far    = Np_GL(4);
% Allocate memory for main matrices
J_mn = zeros(L,M,N,6,6);
G_mn = zeros(L,M,N,3,3);
% Reference cell
r_n = [0.0 , 0.0 , 0.0]';
% Reference distance vector
d = [dx,dy,dz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_far;

[Np,wp,up,vp] = Gauss_2D(N_P);

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Far Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = R_faces;
kko = ko;
ddx = dx;
%

[index_far] = index_assembly(n_medium,L,M,N);
length_index_far = length(index_far);

M_x(:) = index_far(:,1);
M_y(:) = index_far(:,2);
M_z(:) = index_far(:,3);

JJ_mn = zeros(length_index_far,3,3);

parfor i = 1 : length_index_far
    
    mx = M_x(i);
    my = M_y(i);
    mz = M_z(i);
    %
    m = [mx,my,mz];
    r_m = ( (m-1) .* d )';
    
    JJ_mn(i,:,:) = mexCUBATURE_Nop(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
    
end

for i = 1 : length_index_far
    
    mx = M_x(i);
    my = M_y(i);
    mz = M_z(i);
    %
    
    J_mn(mx,my,mz,:,:) = JJ_mn(i,:,:);
    
end
clear JJ_mn;
% parfor mx = 1:L
%     for my = 1:M
%         for mz = 1:N
%             m = [mx,smy,mz];
%             r_m = ( (m-1) .* d )';
%                         
% %             J_mn(mx,my,mz,:,:) = CUBATURE_2D_2D(Np_1D_far,r_m,r_n,RR,kko);
% %             J_mn(mx,my,mz,:,:) = CUBATURE_Nop_sym_far(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
%             J_mn(mx,my,mz,:,:) = mexCUBATURE_Nop(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
%                      
%         end
%     end
% end

Time_far = toc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Medium Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_medium;

[Np,wp,up,vp] = Gauss_2D(N_P);

parfor mx = 1:n_medium
    for my = 1:n_medium
        for mz = 1:n_medium
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
                      
%             J_mn(mx,my,mz,:,:) = CUBATURE_2D_2D(Np_1D_medium,r_m,r_n,RR,kko);
%             J_mn(mx,my,mz,:,:) = CUBATURE_Nop_sym_near(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
            J_mn(mx,my,mz,:,:) = mexCUBATURE_Nop(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);

                     
        end
    end
end

Time_medium = toc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Nearby Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_near;

[Np,wp,up,vp] = Gauss_2D(N_P);
% Compute singular terms via DEMCEM
[I_ST] = DEMCEM_ST(Np_DEMCEM);
[I_EAc,I_EAo] = DEMCEM_EA(Np_DEMCEM);
[I_VAc,I_VAo] = DEMCEM_VA(Np_DEMCEM);
parfor mx = 1:2
    for my = 1:2
        for mz = 1:2
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
            
%             [I_mn_non_corrected] = CUBATURE_2D_2D(Np_1D_near,r_m,r_n,RR,kko);
            [I_mn_non_corrected] = CUBATURE_Nop_sym_near(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);

            
            J_mn(mx,my,mz,:,:) = Nearby_Correct(m,I_mn_non_corrected,I_ST,I_EAc,I_EAo,I_VAc,I_VAo);
 
        end
    end
end

Time_near = toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Evaluate G_mn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
tic

[P_vector] = Vector_Projection();
%
parfor mx = 1:L
    for my = 1:M
        for mz = 1:N
            JJ = squeeze(J_mn(mx,my,mz,:,:));
            for p = 1:3                           
                for q = 1:3
                    PP = squeeze(P_vector(:,:,p,q));
                    G_mn(mx,my,mz,p,q) = sum(sum(JJ .* PP));
                    
                end
            end
                     
        end
    end
end       

Time_Projection = toc

G_mn = ( 1.0 / (4*pi) ) * G_mn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%