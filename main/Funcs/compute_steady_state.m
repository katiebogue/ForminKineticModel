function pi = compute_steady_state(Q)
% COMPUTE_STEADY_STATE  computes the fraction of time spent in each state
% of the transition matrix at steady state
    %
    %   pi=COMPUTE_STEADY_STATE(Q) computes the fraction of time spent in 
    %        each state of the transition matrix, Q, at steady state
    %
    %   Input:
    %       Q    : (sparse matrix) transition matrix, NxN
    %
    %   Output:
    %       pi : matrix of fraction of time spent in each state
    % 
    %   See also FORMIN, FORMINTRANSITIONMAT, PMULTIBIND.

    % Q: Transition rate matrix (NxN matrix)

    % Transpose the matrix Q
    Qt = Q';

    % Find the left eigenvector corresponding to the smallest eigenvalue
    [V, D] = eigs(Qt, 1, 'SM');

    % Extract the eigenvector corresponding to the smallest eigenvalue
    pi = V(:, 1);

    % Normalize the steady-state vector to sum to 1
    pi = pi / sum(pi);

    % Ensure it's a real vector (eliminate any potential complex parts)
    pi = real(pi);
end