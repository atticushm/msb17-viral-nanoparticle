function [nodes,N] = get_virus_nodes(virus_length, n)
%GET_VIRUS_NODES 
%GET_VIRUS_NODES calculates the nodes along a virion constructed from 10
%slender virions.
%   ______________________________________________________________________
%   INPUTS
%   virus_length: contour length of the virion.
%   n: number of desired nodes along the length of a virion.
%   ______________________________________________________________________
%   OUTPUTS
%   nodes: a matrix of nodes describing the surface of the spheroid.
%   N: the number of nodes
%   ______________________________________________________________________

l = virus_length;

int_vir = horzcat(linspace(0,l,n),linspace(-l,0,n));
N = length(int_vir);
int_vir_nodes = vertcat(zeros(1,N),zeros(1,N),int_vir);

nodes1 = get_rotated_nodes(int_vir_nodes,pi/2,0,0);
nodes2 = get_rotated_nodes(int_vir_nodes,0,pi/2,0);
nodes3 = get_rotated_nodes(int_vir_nodes,pi/4,pi/4,0);
nodes4 = get_rotated_nodes(int_vir_nodes,pi/4,-pi/4,0);
nodes5 = get_rotated_nodes(int_vir_nodes,-pi/4,pi/4,0);
nodes6 = get_rotated_nodes(int_vir_nodes,-pi/4,-pi/4,0);

nodes = horzcat(int_vir_nodes,nodes1,nodes2,nodes3,nodes4,nodes5,nodes6);

nodes(:,all(nodes == 0)) = []; % origin (focal point) excluded, should have neglible impact.
N = size(nodes,2);

scatter3(nodes(1,:),nodes(2,:),nodes(3,:))
axis vis3d equal

end
