clear all; clc; close all

% Using Nystrom discretisation.
l = 1;
m = 30; % number of force nodes to put along a single virion.
eps = l/m;
mu = 1;

k = 1;
kbar = 1;

n_v = 1;
v_nodes = vertcat(linspace(0,1,m), zeros(1,m),zeros(1,m));
b0 = v_nodes(1,2)-v_nodes(1,1);
N = length(v_nodes);
figure 
axis vis3d equal
axis([0 1 -1 1 -1 1])
NN = eye(3*N);

h = 0.01;
t = [0:h:1];

x_1 = v_nodes;
scatter3(x_1(1,:),x_1(2,:),x_1(3,:),'ro')

x_curr = [];

for kk = 1:length(t)
    S = get_reg_stokeslets(x_1,x_1,eps,mu);
    F_elas = get_stretch_forces(x_1,b0,k) + get_bend_forces(x_1,kbar);
    
    resist = get_resistance_mat_virion(NN,eps,x_1,N,x_1,N);
    F_hydro = [resist(1,1); resist(2,2); resist(3,3)];
    F = F_elas + kron(F_hydro,ones(1,N));
    F = horzcat(F(1,:),F(2,:),F(3,:))';
    
    x_1 = horzcat(x_1(1,:),x_1(2,:),x_1(3,:))';  
    x_curr = x_1 + h.*S*F;
    x_curr = vertcat(x_curr(1:N)',x_curr(N+1:2*N)',x_curr(2*N+1:end)');
    
    x_1 = x_curr;
    scatter3(x_curr(1,:),x_curr(2,:),x_curr(3,:),'ro')
    axis([0 1 -1 1 -1 1])
    pause(0.01)
end
