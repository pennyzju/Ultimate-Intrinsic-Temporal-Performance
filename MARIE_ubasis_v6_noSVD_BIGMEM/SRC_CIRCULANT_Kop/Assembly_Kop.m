function [IK_mn] = Assembly_Kop(L,M,N,Np_GL,n_medium,ko,dx,dy,dz,R_faces)
%%

% addpath('.\SRC_CIRCULANT_Kop\CUBATURES_Kop')

% Set order of quadrature for each region
% Np_GL = [Np_DEMCEM ; Np_1D_near ; Np_1D_medium ; Np_1D_far];
% Np_DEMCEM    = Np_GL(1);
Np_1D_near   = Np_GL(2);
Np_1D_medium = Np_GL(3);
Np_1D_far    = Np_GL(4);
% Allocate memory for main matrices
IK_mn = zeros(L,M,N,3);
% Reference cell
r_n = [0.0 , 0.0 , 0.0]';
% Reference distance vector
d = [dx,dy,dz];
% modification for "parfor"
RR = R_faces;
kko = ko;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Far Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_far;

[Np,wp,up,vp] = Gauss_2D_Kop(N_P);

parfor mx = 1:L
    for my = 1:M
        for mz = 1:N
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
                        
%             IK_mn(mx,my,mz,:)  = CUBATURE_Kop_far(Np,wp,up,vp,Np,wp,up,vp,r_m,r_n,RR,kko);
            IK_mn(mx,my,mz,:)  = mexCUBATURE_Kop(Np,wp,up,vp,r_m,r_n,RR,kko);
                     
        end
    end
end

Time_far = toc;
fprintf('Time_far         = %d \n',int64(Time_far));
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Medium Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_medium;

[Np,wp,up,vp] = Gauss_2D_Kop(N_P);

parfor mx = 1:n_medium
    for my = 1:n_medium
        for mz = 1:n_medium
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';

%             IK_mn(mx,my,mz,:)  = CUBATURE_Kop_far(Np,wp,up,vp,Np,wp,up,vp,r_m,r_n,RR,kko);
            IK_mn(mx,my,mz,:)  = mexCUBATURE_Kop(Np,wp,up,vp,r_m,r_n,RR,kko);

                     
        end
    end
end

Time_medium = toc;
fprintf('Time_medium      = %d \n',int64(Time_medium));
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Nearby Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_near;

[Np,wp,up,vp] = Gauss_2D_Kop(N_P);

parfor mx = 1:2
    for my = 1:2
        for mz = 1:2
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
                       
            IK_mn(mx,my,mz,:)  = CUBATURE_Kop_near_sym(Np,wp,up,vp,Np,wp,up,vp,r_m,r_n,m,RR,kko);
 
        end
    end
end

Time_near = toc;
fprintf('Time_near        = %d \n',int64(Time_near));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     FINAL OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef = dx^4 / (1i*ko)^2 / (4*pi);
%
IK_mn = coef * IK_mn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%