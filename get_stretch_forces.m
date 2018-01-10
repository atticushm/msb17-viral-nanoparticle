function f = get_stretch_forces (x, b0, k)
%GET_STRETCH_FORCES calculates the stretch force on a string of connected
%nodes
%   _______________________________________________________________________
%   INPUTS:
%   x: nodes presented in a 3xN matrix. Nodes assumed to be equally distributed at equilibrium.
%   b0: equilibrium distance between nodes (scalar value).
%   k: spring constant.
%   _______________________________________________________________________
%   OUTPUTS:
%   f: 3xN matrix of calculated stretching forces.
%   _______________________________________________________________________

bmat = x(:, 2:end) - x(:, 1:end-1);
b0vec = b0*ones(1,length(x(1,:)));

% make norm vector
bb = bmat.*bmat;
b_norm_vec = sqrt(bb(1,:)+bb(2,:)+bb(3,:));

% matrix of b vectors
b_n_vecs = [bmat, zeros(3,1)];
b_nneg1_vecs = [zeros(3,1), bmat];

% vector of b norms
b_n = [b_norm_vec, 1];
b_nneg1 = [1,b_norm_vec];

% coefficients:
b_n_coeff = (b_n-b0vec)./b_n;
b_nneg1_coeff = (b_nneg1-b0vec)./b_nneg1;
b_n_coeff_mat = [b_n_coeff; b_n_coeff; b_n_coeff];
b_nneg1_coeff_mat = [b_nneg1_coeff; b_nneg1_coeff; b_nneg1_coeff];

f = k*b_n_coeff_mat.*b_n_vecs - k*b_nneg1_coeff_mat.*b_nneg1_vecs;

end
 
