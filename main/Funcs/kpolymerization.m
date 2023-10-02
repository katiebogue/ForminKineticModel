function output = kpolymerization(type,kcap,kdel,rcap,rdel,krel)
%KPOLYMERIZATION Computes the rate of elongation for a single PRM
    %   kpoly=
    %   KPOLYMERIZATION(type,kcap,kdel,rcap,rdel,krel) computes the rate of
    %   elongation
    %
    %   Inputs:
    %         type : (string) the elongation rate equation type; either
    %                 "capture", "3st", or "4st"
    %         kcap : (double or FilType) rate of capture for the PRM
    %         kdel : (double or FilType) rate of delivery for the PRM
    %         rcap : (double or FilType) reverse rate of capture for the 
    %                PRM
    %         rdel : (double or FilType) reverse rate of delivery for the 
    %                 PRM
    %         krel : (double or FilType) rate of release for the PRM
    %
    %   Output is either double or FilType (depending on inputs)
    %
    %   See also KPOLY, PRM, FORMIN, FILTYPE.
arguments
    type (1,1) string {mustBeMember(type,{'capture','3st','4st'})}
    kcap {mustBeA(kcap,["double","FilType"])}
    kdel {mustBeA(kdel,["double","FilType"])}
    rcap {mustBeA(rcap,["double","FilType"])}
    rdel {mustBeA(rdel,["double","FilType"])}
    krel {mustBeA(krel,["double","FilType"])}
end
    if type=="capture"
        output=kcap;
    elseif type=="3st"
        output=1/((1/kdel) + ((kdel + rcap)/(kdel*kcap)));
    elseif type=="4st"
        output=1/((1/krel) + ((rdel + krel)/(kdel * krel)) + (((rcap * rdel) + (rcap * krel) + (kdel * krel))/(kcap * kdel * krel)));
    end
end