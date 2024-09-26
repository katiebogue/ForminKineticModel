function frac = pmultibind(formin1)
% PMULTIBIND  computes the fraction of time spent with more than
% one PRM bound for the input formin assuming the 3 state model
    %
    %   frac=PMULTIBIND(formin1) computes fraction of time with multiple simultaneous binding for formin1
    %
    %   Input:
    %       formin1    : Formin object, with rates caluclated according to
    %       the formin's options property
    %
    %   Output:
    %       frac : (FilType) fraction of time spent with more than one PRM 
    %               bound for the single, dimer, and double models of the 
    %               3 state model
    %
    %   Calls formintransitionmat and compute_steady_state
    % 
    %   See also FORMIN, FORMINTRANSITIONMAT, COMPUTE_STEADY_STATE.
    arguments
        formin1 Formin
    end
    % compute the fraction of time the formin has multiple PRMs bound to
    % profilin-actin at once
    
    [Qstruct,nbound]=formintransitionmat(formin1); % transition matrix and number of PRMs bound per state

    frac=FilType;

    frac.single=computefrac(Qstruct.single,nbound.single);
    frac.double=computefrac(Qstruct.double,nbound.double);
    frac.dimer=computefrac(Qstruct.dimer,nbound.dimer);
end

function frac = computefrac(Q,nbound)
    pi = compute_steady_state(Q); % fraction of time spend in each state at steady state
    
    index=nbound>1;
    
    frac=sum(pi(index));
end