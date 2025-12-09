function pval=pr(n1,fh1length,FH2size,k,x,y,type,scaledTF)
% PR calculate the probability density of a PRM at a specific location
% using equations
    %   pval=
    %   PR(n1,fh1length,FH2size,k,x,y,type,scaledTF) calculates the probability
    %   density at the x and y coordinates of the PRM at location n1 on a FH1 of length fh1length,
    %   with FH2 size FH2 size, for the specified filament type
    %
    %   Inputs:
    %         n1        : location of PRM, number of amino acids from the FH2
    %         fh1length : FH1 length
    %         FH2size   : distance between the C-terminal amino acids
    %                         of the FH1
    %         k         : kuhn lengths to use (should set to 1)
    %         x         : x coordinate of delivery location (the FH2 plane)
    %         y         : y coordinate of delivery location (perpendicular to the FH2)
    %         type      : (string) "dimer," "double," or "ratio." Note that
    %                     the double value is the same as single, since they use the
    %                     same equation.
    %         scaledTF  : (bool) whether or not to scale the results to be
    %         in uM-1 rather than kuhn length units (default is false)
    %
    %   Outputs:
    %         pval : (double) the probability density
    %   
    %   See also PRM, FORMIN.

    arguments
        n1
        fh1length
        FH2size
        k
        x
        y
        type
        scaledTF= false
    end
    nN=n1+2*(fh1length-n1);
    nprime=((1./n1)+(1./nN)).^-1;
    b=FH2size;
    x2=b;
    y2=0;
    dotprod=(x.*x2)+(y.*y2);
    pval_dimer= (...
        (3./(2*pi.*nprime.*(k.^2))).^(3/2) .* ...
        exp( (3.*(b^2)) ./ (2*(k^2).*(n1+nN)) ).*...
        exp( (-3./(2.*(k^2))) .* ( ((x.^2+y.^2)./n1) + ((x.^2+y.^2+b^2-2.*dotprod)./nN) ) )...
    );
    pval_double=( ( 3./(2*pi.*n1.*(k^2)) ).^(3/2)) .* exp( -3.*(x.^2+y.^2)./(2.*n1.*(k^2)) ) ;

    if scaledTF
        pval_dimer=pval_dimer.*6.1503*10^7;
        pval_double=pval_double.*6.1503*10^7;
    end

    if type=="dimer"
        pval=pval_dimer;
    elseif type=="double"
        pval=pval_double;
    elseif type=="ratio"
        pval=pval_dimer./pval_double;
    end
end