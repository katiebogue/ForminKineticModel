function output = kpolymerization(type,kcap,kdel,rcap,rdel,krel)
%KPOLYMERIZATION Computes the rate of elongation for a single PRM
    %   kpoly=
    %   KPOLYMERIZATION(type,kcap,kdel,rcap,rdel,krel) computes the rate of
    %   elongation
    %
    %   Inputs:
    %         type : (string) the elongation rate equation type; either
    %                 "capture", "3st", or "4st"
    %         kcap : (double or Submodel) rate of capture for the PRM
    %         kdel : (double or Submodel) rate of delivery for the PRM
    %         rcap : (double or Submodel) reverse rate of capture for the 
    %                PRM
    %         rdel : (double or Submodel) reverse rate of delivery for the 
    %                 PRM
    %         krel : (double or Submodel) rate of release for the PRM
    %
    %   Output is either double or submodel (depending on inputs)
    %
    %   See also KPOLY, PRM, FORMIN, SUBMODEL.
arguments
    type (1,1) string {mustBeMember(type,{"capture","3st","4st"})}
    kcap {mustBeA(kcap,["double","Submodel"])}
    kdel {mustBeA(kdel,["double","Submodel"])}
    rcap {mustBeA(rcap,["double","Submodel"])}
    rdel {mustBeA(rdel,["double","Submodel"])}
    krel {mustBeA(krel,["double","Submodel"])}
end
    if type=="capture"
        output=kcap;
    elseif type=="3st"
        output=1/((1/kdel) + ((kdel + rcap)/(kdel*kcap)));
    elseif type=="4st"
        output=1/((1/krel) + ((rdel + krel)/(kdel * krel)) + (((rcap * r_paf_rev) + (rcap * krel) + (kdel * rdel))/(kcap * kdel * krel)));
    end
end