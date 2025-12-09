function [c_PA,c_actin_free,c_profilin_free]=calcCPA(c_actin,c_profilin,kd)
% CALCCPA calculate the [profilin-actin] from the initial [actin], initial
% [profilin], and the kd
    %   [c_PA,c_actin_free,c_profilin_free]= CALCCPA(c_actin,c_profilin,kd) 
    %
    %   Inputs:
    %         c_actin    : (double) the initial [actin]
    %         c_profilin : (double) the initial [profilin]
    %         kd         : (double) dissociation constant
    %
    %   Ensure units are the same for all 3 inputs; Output units are same
    %   as input
    %
    %   Outputs:
    %         c_PA            : (double) the [profilin-actin]
    %         c_actin_free    : (double) the free [actin]
    %         c_profilin_free : (double) the free [profilin]
    %   
    %   See also PRM, FORMIN.

% kd = (c_actin_free .* c_profilin_free)./c_PA

c_PA=((c_actin+c_profilin+kd)-sqrt((c_actin+c_profilin+kd).^2-4.*(c_actin.*c_profilin)))./2;

c_actin_free=c_actin-c_PA;

c_profilin_free=c_profilin-c_PA;
end


