function frac = pmultibind(formin1)
    arguments
        formin1 Formin
    end
    % compute the fraction of time the formin has multiple PRMs bound to
    % profilin-actin at once
    
    [Qstruct,nbound]=formintransitionmat(formin1);

    frac=FilType;

    frac.single=computefrac(Qstruct.single,nbound.single);
    frac.double=computefrac(Qstruct.double,nbound.double);
    frac.dimer=computefrac(Qstruct.dimer,nbound.dimer);
end

function frac = computefrac(Q,nbound)
    pi = compute_steady_state(Q);
    
    index=nbound>1;
    
    frac=sum(pi(index));
end