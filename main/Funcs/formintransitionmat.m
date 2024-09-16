function [Qstruct,nbound]=formintransitionmat(formin1)
    arguments
        formin1 Formin
    end
    Qstruct=FilType;
    nbound=FilType;
    N=formin1.PRMCount;
    kons=FilType(zeros(1,N),zeros(1,2*N),zeros(1,2*N));
    koffs=kons;

    for k=1:N
        kcap=formin1.PRMList(k).kcap;
        rcap=formin1.PRMList(k).rcap;
        kdel=formin1.PRMList(k).kdel;
        
        kons.single(k)=kcap.single;
        kons.double(k)=kcap.double.a;
        kons.double(k+N)=kcap.double.b;
        kons.dimer(k)=kcap.dimer.a;
        kons.dimer(k+N)=kcap.dimer.b;

        % Use sum of off-rates (parallel pathways assumption)
        koffs.single(k)=kdel.single+rcap;
        koffs.double(k)=kdel.double.a+rcap;
        koffs.double(k+N)=kdel.double.b+rcap;
        koffs.dimer(k)=kdel.dimer.a+rcap;
        koffs.dimer(k+N)=kdel.dimer.b+rcap;
    end
    [Qstruct.single,nbound.single]=generate_transition_matrix(N,kons.single,koffs.single);
    [Qstruct.double,nbound.double]=generate_transition_matrix(2*N,kons.double,koffs.double);
    [Qstruct.dimer,nbound.dimer]=generate_transition_matrix(2*N,kons.dimer,koffs.dimer);

    function [Q,nbound] = generate_transition_matrix(N, k_on, k_off)
        % N: Number of binding sites
        % k_on: Vector of on-rates [k_on1, k_on2, ..., k_onN]
        % k_off: Vector of off-rates [k_off1, k_off2, ..., k_offN]
    
        % Number of states (2^N, since each site can be bound or unbound)
        num_states = 2^N;
    
        % Generate the states as binary vectors
        states = dec2bin(0:num_states-1) - '0'; % Each row is a state
        nbound=sum(states,2);
    
        % Preallocate arrays for row indices, column indices, and values
        max_nonzeros = num_states * (N + 1);
        row_indices = zeros(max_nonzeros, 1);
        col_indices = zeros(max_nonzeros, 1);
        values = zeros(max_nonzeros, 1);
        
        % Initialize a counter to keep track of the position in the preallocated arrays
        idx = 1;
    
        % Loop over all states
        for i = 1:num_states
            for j = 1:N
                if states(i, j) == 0
                    % If the binding site j is unbound in state i, bind it
                    new_state = states(i, :);
                    new_state(j) = 1;
                    new_state_index = bin2dec(num2str(new_state)) + 1;
    
                    % Record the transition rate
                    row_indices(idx) = i;
                    col_indices(idx) = new_state_index;
                    values(idx) = k_on(j);
                    idx = idx + 1;
                else
                    % If the binding site j is bound in state i, unbind it
                    new_state = states(i, :);
                    new_state(j) = 0;
                    new_state_index = bin2dec(num2str(new_state)) + 1;
    
                    % Use the sum of the off-rates
                    row_indices(idx) = i;
                    col_indices(idx) = new_state_index;
                    values(idx) = k_off(j);
                    idx = idx + 1;
                end
            end
    
            % Set the diagonal element (rate of leaving state i)
            row_indices(idx) = i;
            col_indices(idx) = i;
            values(idx) = -sum(values(idx-N:idx-1));
            idx = idx + 1;
        end
    
        % Remove any unused preallocated space
        row_indices = row_indices(1:idx-1);
        col_indices = col_indices(1:idx-1);
        values = values(1:idx-1);
    
        % Create the sparse matrix
        Q = sparse(row_indices, col_indices, values, num_states, num_states);
    end
end