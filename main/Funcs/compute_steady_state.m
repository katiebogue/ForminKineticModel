function pi = compute_steady_state(Q)
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