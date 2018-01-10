function [G_resist_calc, hf, hq] = get_resistance_mat_virion(nn_mat,eps,f_nodes,N,q_nodes,Q)
%GET_RESISTANCE_MAT computes the approximate grand resistance matrices for a sphere in stokes flow.
%   ______________________________________________________________________
%   This function uses the functions:
%       get_ellipsoidal_nodes
%       get_reg_stokeslets
%       get_rotated_nodes
%       get_nn_mat
%       get_grid_spacing
% 
%   The problem solved is dimensionless (ie mu=1). Nystrom discretisation
%   (ie equal force and quadrature discretisations) is used.
%   ______________________________________________________________________
%   INPUTS
%   n_v: number of virions joined together to form the nanoparticle.
%   eps: regularisation parameter epsilon
%   f_dim: vector of dimensions of cube for force nodes projection.
%   q_dim: vector of dimensions of cube for quadrature nodes projection.
%   ______________________________________________________________________
%   OUTPUTS
%   G_resist_calc: the Grand Resistance Matrix for the given object.
%   hf: max force discretisation spacing.
%   hq: max quadrature discretisation spacing.
%   ______________________________________________________________________

mu = 1;              

prescribed = -0.01.*eye(6);            % unit translations and unit angular velocities
resist_T = zeros(3);
resist_R = zeros(3);
resist_P1 = zeros(3);
resist_P2 = zeros(3);

[hf, hq] = get_grid_spacing(f_nodes,q_nodes);

% Calculate stokeslets and thus reduced LHS using NN discretisation.
S = get_reg_stokeslets(q_nodes,f_nodes,eps,mu);
A = S*nn_mat;

clear kk
for kk = 3

    u_int = prescribed(kk,:);

    U_trans = horzcat(u_int(1).*ones(1,N),u_int(2).*ones(1,N),u_int(3).*ones(1,N));
    Omg = vertcat(u_int(4).*ones(1,N),u_int(5).*ones(1,N),u_int(6).*ones(1,N));  
    U_rot = cross(Omg,f_nodes,1);
    U_rot = horzcat(U_rot(1,:),U_rot(2,:),U_rot(3,:));
    U = U_trans + U_rot;

%   Forces.
    b = A\U';
    for n = 1:3*N
        f(n) = sum(dot(b(n).*ones(1,3*Q),nn_mat(:,n)));
    end
    F_sums = zeros(1,3);
    F_sums(1) = sum(f(1:N));
    F_sums(2) = sum(f(N+1:2*N));
    F_sums(3) = sum(f(2*N+1:3*N));

%   Moments.
    f = reshape(f,[N,3])';
    M = cross(f_nodes,f);   
    M_sums = zeros(1,3);
    M_sums(1) = sum(M(1,:));
    M_sums(2) = sum(M(2,:));
    M_sums(3) = sum(M(3,:));
    
    if max(kk == [1 2 3]) == 1
        resist_T(:,kk) = F_sums;
        resist_P1(:,kk) = M_sums;
    else
        resist_R(:,kk-3) = M_sums;
        resist_P2(:,kk-3) = F_sums;
    end
    
end

G_resist_calc = vertcat(horzcat(resist_T, resist_P1),horzcat(resist_P2, resist_R));

end
