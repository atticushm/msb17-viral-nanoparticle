clear all; clc; close all

% Using Nystrom discretisation.
l = 1;
m = 15; % number of force nodes to put along a single virion.
eps = l/m;
mu = 1;

k = 1;
kbar = 1;

% Location of end points of N virions is given by Euan in "virion_endpoints.txt". We use this to
% construct the nodes for both force and quadrature, as well as nearest
% neighbour matrices all at once to ensure force nodes only connect the
% quadrature nodes from the same virion.
endpoints = dlmread('virion_endpoints_n26.txt')';
n_v = size(endpoints,2);
nodes = [];
F_mat = [];

for kk = 1:n_v
    fx = linspace(0,endpoints(1,kk),m);
    fy = linspace(0,endpoints(2,kk),m);
    fz = linspace(0,endpoints(3,kk),m);
    
    fx = 0.5*(fx(1:end-1) + fx(2:end));
    fy = 0.5*(fy(1:end-1) + fy(2:end));
    fz = 0.5*(fz(1:end-1) + fz(2:end));
       
    f_n = vertcat(fx,fy,fz);   
    f_n(:,all(f_n == 0)) = [];
    M(kk) = length(f_n);
    
    b0 = f_n(end-3) - f_n(end-4);
    F1 = get_bend_forces(f_n,kbar);
    F2 = get_stretch_forces(f_n,b0,k);
    F_n = F1+F2;
    F_n(:,1) = zeros(3,1);
    
    F_mat = horzcat(F_mat,F_n);
    nodes = horzcat(nodes, f_n);
end
N = size(nodes,2);
h = 0.01;
t = [0:h:1];
NN = eye(3*N);

x_2 = horzcat(nodes(1,:),nodes(2,:),nodes(3,:));
scatter3(x_2(1:N),x_2(N+1:2*N),x_2(2*N+1:end),'.')

G(length(t)) = struct('cdata',[],'colormap',[]);
figure
axis vis3d equal
axis([-2 2 -2 2 -8 2])
xlabel('x')
ylabel('y')
zlabel('z')
title('t = 0')
drawnow
G(1) = getframe(gcf);
pause(0.01)

% Forward Euler method for first time step.
x_2 = vertcat(x_2(1:N),x_2(N+1:2*N),x_2(2*N+1:end));
S = get_reg_stokeslets(x_2,x_2,eps,mu);

resist1 = get_resistance_mat_virion(NN,eps,x_2,N,x_2,N);
F_hydro = [resist1(1,1); resist1(2,2); resist1(3,3)];
F_hydro = kron(F_hydro,ones(1,N));
F_hydro = horzcat(F_hydro(1,:),F_hydro(2,:),F_hydro(3,:));

F_vec = horzcat(F_mat(1,:),F_mat(2,:),F_mat(3,:))';
F = F_vec + F_hydro';
        
x_2 = horzcat(nodes(1,:),nodes(2,:),nodes(3,:))';
x_1 = x_2 + h.*S*F;
scatter3(x_1(1:N),x_1(N+1:2*N),x_1(2*N+1:end),'.')
axis vis3d equal
axis([-2 2 -2 2 -8 2])
drawnow
G(2) = getframe(gcf);
pause(0.01)

% Adam Bashford on remaining time steps.
for n = 3:length(t)
    x_1 = vertcat(x_1(1:N)',x_1(N+1:2*N)',x_1(2*N+1:end)');
    x_2 = vertcat(x_2(1:N)',x_2(N+1:2*N)',x_2(2*N+1:end)');
    F_1 = [];
    F_2 = [];
    for kk = 1:n_v
        nodes_1 = x_1(:,((kk-1)*M(kk)+1):kk*M(kk));
        nodes_2 = x_2(:,((kk-1)*M(kk)+1):kk*M(kk));
        
        F1 = get_bend_forces(nodes_1,kbar) + get_stretch_forces(nodes_1,b0,k);
        F2 = get_bend_forces(nodes_2,kbar) + get_stretch_forces(nodes_2,b0,k);
        
        F_1 = horzcat(F_1,F1);
        F_2 = horzcat(F_2,F2);
    end   
    
    resist1 = get_resistance_mat_virion(NN,eps,x_1,N,x_1,N);
    F_hydro1 = [resist1(1,1); resist1(2,2); resist1(3,3)];
    resist2 = get_resistance_mat_virion(NN,eps,x_2,N,x_2,N);
    F_hydro2 = [resist2(1,1); resist2(2,2); resist2(3,3)];
    
    F_1 = F_1 + kron(F_hydro1,ones(1,N));
    F_2 = F_2 + kron(F_hydro2,ones(1,N));
    
    F_1 = horzcat(F_1(1,:),F_1(2,:),F_1(3,:))';
    F_2 = horzcat(F_2(1,:),F_2(2,:),F_2(3,:))';    
    
    Fun1 = get_reg_stokeslets(x_1,x_1,eps,mu)*F_1;
    Fun2 = get_reg_stokeslets(x_2,x_2,eps,mu)*F_2;
    
    x_1 = horzcat(x_1(1,:),x_1(2,:),x_1(3,:))';
    x_2 = horzcat(x_2(1,:),x_2(2,:),x_2(3,:))';  
    x_curr = x_1 + (h/2).*(3.*Fun1 - Fun2);
    
    x_1 = x_curr;
    x_2 = x_1;
    
    scatter3(x_curr(1:N),x_curr(N+1:2*N),x_curr(2*N+1:end),'.')
    axis vis3d equal
    axis([-2 2 -2 2 -8 2])
    str = sprintf('t = %g',t(n));
    title(str)
    drawnow
    G(n) = getframe(gcf);
    pause(0.01)
end

video = VideoWriter('Nanoparticle n16 sedimenting.avi','Motion JPEG AVI');
video.FrameRate = 10;
open(video)
writeVideo(video,G)
close(video)

