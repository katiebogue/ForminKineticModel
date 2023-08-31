function output = kcapture(kcap_scale,p_occ,c_PA)
%KCAPTURE Computes the rate of capture for a single PRM
    %   kcap=
    %   KCAPTURE(kcap_scale,p_occ,c_PA) computes the rate of capture
    %
    %   Inputs:
    %         kcap_scale : the capture rate constant (should already be
    %                      scaled based on additonal parameters)
    %         p_occ      : the occlusion probability of the PRM
    %         c_PA       : concentration of profilin-actin   
    %
    %   See also KPOLY, PRM, FORMIN.
output=kcap_scale*c_PA*(1-p_occ);
end