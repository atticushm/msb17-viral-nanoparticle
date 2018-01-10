function nn_virion = get_nn_virion(coarse_nodes, fine_nodes)
%   ______________________________________________________________________
%   INPUTS
%   coarse_nodes: list of nodes in the course discretisation.
%   fine_nodes: list of nodes in the finer discretisation.
%
%   For example, a list of 3-dimensional nodes should be inputted as a 3xN
%   matrix, with the first, second, and third rows containing the x, y, and
%   z values of each node respectively.
%   ______________________________________________________________________
%   OUTPUTS
%   nn_virion: 3Qx3N nearest-neighbour matrix
%   ______________________________________________________________________

x = coarse_nodes;
X = fine_nodes;

dim = size(x,1);
N = size(x,2);
Q = size(X,2);

diff_mat_x = repmat(x(1,:),Q,1) - repmat(X(1,:)',1,N);
diff_mat_y = repmat(x(2,:),Q,1) - repmat(X(2,:)',1,N);
diff_mat_z = repmat(x(3,:),Q,1) - repmat(X(3,:)',1,N);

diff_mat = sqrt(diff_mat_x.^2 + diff_mat_y.^2 + diff_mat_z.^2);
[~, nhat] = min(diff_mat,[],2);

nn_virion = zeros(Q,N);
for q = 1:Q
    nn_virion(q,nhat(q)) = 1;
end

nn_virion = kron(speye(dim),nn_virion);

end

