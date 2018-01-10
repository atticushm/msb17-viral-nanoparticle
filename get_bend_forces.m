function f = get_bend_forces(x, kbar)
%GET_BEND_FORCES calculates the bending force on a string of connected
%nodes
%   _______________________________________________________________________
%   INPUTS:
%   x: nodes presented in a 3xN matrix. Nodes assumed to be equally distributed at equilibrium.
%   kbar: bending constant.
%   _______________________________________________________________________
%   OUTPUTS:
%   f: 3xN matrix of calculated stretching forces.
%   _______________________________________________________________________

bmat = x(:, 2:end) - x(:, 1:end-1);

b_n = [bmat, zeros(3,1)];
b_nneg1 = [zeros(3,1), bmat];
b_nneg2 = [zeros(3,1), b_nneg1(:, 1:end-1)];
b_npl1 = [b_n(:, 2:end), zeros(3,1)];

f = 3*b_n - 3*b_nneg1 + b_nneg2 - b_npl1;
f = kbar*f;

end
    
    
