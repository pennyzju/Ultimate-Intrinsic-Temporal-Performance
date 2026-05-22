function [IK_mn] = CUBATURE_Kop_near(Np,wp,up,vp,Nq,wq,uq,vq,r_m,r_n,m)
%%
global R_faces
global ko
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                1D cubature's number of points                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_mn = zeros(6,6,3);
J_mn = zeros(6,6);
%
Rm   = repmat(r_m,1,Np);
Rn   = repmat(r_n,1,Nq);
%
Rmn_x = zeros(Np,Np);
Rmn_y = zeros(Np,Np);
Rmn_z = zeros(Np,Np);
%%
for i_obs = 1:6
    
    rp1 = R_faces(:,1,i_obs);
    rp2 = R_faces(:,2,i_obs);
    % rp3 = R_faces(:,3,i_obs);
    rp4 = R_faces(:,4,i_obs);
    
    for i_src = 1:6
        
        rq1 = R_faces(:,1,i_src);
        rq2 = R_faces(:,2,i_src);
        % rq3 = R_faces(:,3,i_src);
        rq4 = R_faces(:,4,i_src);
        %                    Opservation face P                            
        a_p = rp2-rp1;
        b_p = rp4-rp1;
        %                         Source face P                            
        a_q = rq2-rq1;
        b_q = rq4-rq1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 4D Numerical Cubature                           %
        %                    Fully vectorized                             % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Weights
        W = wp*wq';
        % Observation points
        Rp_1 = repmat(rp1,1,Np);
        rrp = Rp_1 + a_p * up' + b_p * vp' + Rm;
        % Source points
        Rq_1 = repmat(rq1,1,Nq);
        rrq = Rq_1 + a_q * uq' + b_q * vq' + Rn;
        %
        
        for ii=1:Np
            Rmn_x(ii,:) = rrp(1,ii)-rrq(1,:);
            Rmn_y(ii,:) = rrp(2,ii)-rrq(2,:);
            Rmn_z(ii,:) = rrp(3,ii)-rrq(3,:);
        end
        %
        Rpq = sqrt(Rmn_x.^2 + Rmn_y.^2 + Rmn_z.^2);
        % Singular kernel      
        G = (1-exp(-1i*ko*Rpq) .* (1 + 1i*ko*Rpq) )./ Rpq.^3;
        % Integral
        I_mn(i_obs,i_src,1) = sum(sum(W .* G .* Rmn_x));
        I_mn(i_obs,i_src,2) = sum(sum(W .* G .* Rmn_y));
        I_mn(i_obs,i_src,3) = sum(sum(W .* G .* Rmn_z));
    end
end
% Coincident surfaces' contribution = 0
[I_mn] = Nearby_Correct_Kop(m,I_mn);
%   (uninormal_k,I_mn)
% uninormal_k = [-x +x -y +y -z +z]';
J_mn(1,:) = -I_mn(1,:,1);
J_mn(2,:) =  I_mn(2,:,1);
J_mn(3,:) = -I_mn(3,:,2);
J_mn(4,:) =  I_mn(4,:,2);
J_mn(5,:) = -I_mn(5,:,3);
J_mn(6,:) =  I_mn(6,:,3);
%
IK_mn_x = sum(J_mn(:,2) - J_mn(:,1));
IK_mn_y = sum(J_mn(:,4) - J_mn(:,3));
IK_mn_z = sum(J_mn(:,6) - J_mn(:,5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FINAL OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IK_mn = [IK_mn_x IK_mn_y IK_mn_z]'; % without dx^2 = area of panel = Jacobian