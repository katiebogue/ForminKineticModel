function poccval=poccc(n1,fh1length,type, k,R,epsilon)
% POCCC calculate the occlusion probability of a PRM at a specific location
% using equations, only works for single or double filament with FH2 size=0
    %
    %   Inputs:
    %         n1        : location of PRM, number of amino acids from the FH2
    %         fh1length : FH1 length
    %         type      : (string) "single" or "double"
    %         k         : kuhn lengths to use (default is 1)
    %         R         : radius of sphere (default is 2.25 kuhn lengths)
    %         epsilon   : starting position (default is 1)
    %
    %   Outputs:
    %         poccval : (double) the occlusion probability
    %   
    %   See also PRM, FORMIN.

    arguments
        n1
        fh1length
        type
        k=1
        R=2.25
        epsilon=1
    end
    
    if type=="single"
        poccval= 1- ( ...
            ( 1-((R/(R+epsilon))*erfc((epsilon/k)*(3/(fh1length-n1))^0.5)) )* ... // occlusion from PRM to NTD
            ( 1-((R/(R+epsilon))*erfc((epsilon/k)*(3/n1)^0.5)) ) ... // occlusion from PRM to FH2
            );
    elseif type=="double"
        poccval= 1- ( ...
            ( 1-((R/(R+epsilon))*erfc((epsilon/k)*(3/(fh1length-n1))^0.5)) )* ... // occlusion from PRM to cis NTD
            ( 1-((R/(R+epsilon))*erfc(((epsilon/k)*3/(n1+fh1length))^0.5)) ) ... // occlusion from PRM to trans NTD
            );
    end

end