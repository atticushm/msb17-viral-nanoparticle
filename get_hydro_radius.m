function a = get_hydro_radius(filename)

l = 1;
q = 150; % number of stokeslets to put along a single virion.
n = 10; % number of force nodes to put along a single virion.
eps = l/q;

% Location of end points of N virions is given by Euan in "virion_endpoints.txt". We use this to
% construct the nodes for both force and quadrature, as well as nearest
% neighbour matrices all at once to ensure force nodes only connect the
% quadrature nodes from the same virion.
endpoints = dlmread(filename)';
n_v = size(endpoints,2);
f_nodes = [];
q_nodes = [];
NN = [];

for kk = 1:n_v
    fx = linspace(0,endpoints(1,kk),n);
    fy = linspace(0,endpoints(2,kk),n);
    fz = linspace(0,endpoints(3,kk),n);
    qx = linspace(0,endpoints(1,kk),q);
    qy = linspace(0,endpoints(2,kk),q);
    qz = linspace(0,endpoints(3,kk),q);
    
    fx = 0.5*(fx(1:end-1) + fx(2:end));
    fy = 0.5*(fy(1:end-1) + fy(2:end));
    fz = 0.5*(fz(1:end-1) + fz(2:end));
       
    f_n = vertcat(fx,fy,fz);
    q_n = vertcat(qx,qy,qz);
    
    f_n(:,all(f_n == 0)) = [];
    q_n(:,all(q_n == 0)) = [];
    
    nn = get_nn_virion(f_n,q_n);
    NN = blkdiag(NN, nn);
    f_nodes = horzcat(f_nodes, f_n);
    q_nodes = horzcat(q_nodes, q_n); 
end

% N and Q are the number of force and quadrature nodes respectively.
N = size(f_nodes,2);
Q = size(q_nodes,2);

% Calculate grand resistance matrix by method of regularised stokeslets.
resist_mat = get_resistance_mat_virion(NN,eps,f_nodes,N,q_nodes,Q);

F = resist_mat(1:3,1:3);
T = resist_mat(4:6,4:6);

a = mean(diag((1/(6*pi)).*F));

% aT = mean(diag(((1/(8*pi)).*T).^(1/3)));

hold on
scatter3(q_nodes(1,:),q_nodes(2,:),q_nodes(3,:))
scatter3(f_nodes(1,:),f_nodes(2,:),f_nodes(3,:),'filled')
axis vis3d equal
hold off

end
