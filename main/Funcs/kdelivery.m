function output = kdelivery(kdel_scale,p_occ0,p_r)
%KDELIVERY Computes the rate of delivery for a single PRM
    %   kdel=
    %   KDELIVERY(kcap_scale,p_occ0,c_PA) computes the rate of capture
    %
    %   Inputs:
    %         kdel_scale : the delivery rate constant (should already be
    %                      scaled based on additonal parameters)
    %         p_occ0     : the occlusion probability of the barbed end
    %                      (delivery site)
    %         p_r        : probability density of the PRM at the barbed 
    %                      end (delivery site)  
    %
    %   See also KPOLY, PRM, FORMIN.
output=kdel_scale*(1-p_occ0)*(1.0e33*p_r/(27*6.022e23));
end