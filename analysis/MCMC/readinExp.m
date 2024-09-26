function out_struct=readinExp(Exp)
%READINEXP helper function for MCMCParamfit, reads information from
%Experiment object into a useful struct
%
%   out_struct = READINEXP(Exp) create output structure from values in Exp,
%   to be used in MCMCParamfit
%
%   Inputs:
%       Exp : Experiment containing the fit data and options to be used in
%       the MCMC fit
%       
%   Outputs:
%       out_struct : structure with the following fields: 
%           data - Exp.data 
%           rates - struct with rate values using 1 for all rate constants
%           opts - Exp.opts
%           fh1sizes - cell array of FH1 sizes of formins in data (indexed
%               by each PRM)
%           prmlocs - cell array of PRM locations of formins in data (indexed
%               by each PRM)
%   
% See also MCMCPARAMFIT, EXPERIMENT, OPTIONS.
    % read in formin info 
    ogtable=struct('k_cap',Exp.opts.k_cap,'k_del',Exp.opts.k_del,'r_cap',Exp.opts.r_cap,'r_cap_exp',Exp.opts.r_cap_exp,'r_del',Exp.opts.r_del,'k_rel',Exp.opts.k_rel);
    Exp.opts.k_cap= 1;
    Exp.opts.k_del=1; 
    Exp.opts.r_cap=1; 
    Exp.opts.r_del=1; 
    Exp.opts.k_rel=1;
    Exp.opts.r_cap_exp=1;

    data=Exp.data;
    L=length(data);

    k_capbase=cell([L 1]);
    k_delbase=cell([L 1]);
    r_capbase=cell([L 1]);
    r_delbase=cell([L 1]);
    k_relbase=cell([L 1]);

    fh1sizes=cell([L 1]);
    prmlocs=cell([L 1]);

    for i=1:length(data)
        formin=data(i).formin;
        kcaps=zeros(formin.PRMCount,5);
        kdels=zeros(formin.PRMCount,5);
        rcaps=zeros(formin.PRMCount,5);
        rdels=zeros(formin.PRMCount,5);
        krels=zeros(formin.PRMCount,5);

        fh1s=zeros(formin.PRMCount,1);
        prms=zeros(formin.PRMCount,1);
        for j=1:formin.PRMCount
            prm1=formin.PRMList(1,j);

            fh1s(j,:)=prm1.fh1length;
            prms(j,:)=prm1.dist_FH2;

            if class(prm1.kcap)=="FilType"
                kcaps(j,:)=prm1.kcap.fil2array;
            else
                kcaps(j,:)=prm1.kcap;
            end
            if class(prm1.kdel)=="FilType"
                kdels(j,:)=prm1.kdel.fil2array;
            else
                kdels(j,:)=prm1.kdel;
            end
            if class(prm1.rcap)=="FilType"
                rcaps(j,:)=prm1.rcap.fil2array;
            else
                rcaps(j,:)=prm1.rcap;
            end
            if class(prm1.rdel)=="FilType"
                rdels(j,:)=prm1.rdel.fil2array;
            else
                rdels(j,:)=prm1.rdel;
            end
            if class(prm1.krel)=="FilType"
                krels(j,:)=prm1.krel.fil2array;
            else
                krels(j,:)=prm1.krel;
            end
        end
        k_capbase{i,1}=kcaps;
        k_delbase{i,1}=kdels;
        r_capbase{i,1}=rcaps;
        r_delbase{i,1}=rdels;
        k_relbase{i,1}=krels;

        fh1sizes{i,1}=fh1s;
        prmlocs{i,1}=prms;
    end
    % get values with all parameters set to 1

    % Return to original opts values 
    Exp.opts.k_cap=ogtable.k_cap;
    Exp.opts.k_del=ogtable.k_del;
    Exp.opts.r_cap=ogtable.r_cap;
    Exp.opts.r_cap_exp=ogtable.r_cap_exp;
    Exp.opts.r_del=ogtable.r_del;
    Exp.opts.k_rel=ogtable.k_rel;

    out_struct.data=data;
    out_struct.data = rmfield(out_struct.data, 'formin');

    rates.k_capbase=k_capbase;
    rates.k_delbase=k_delbase;
    rates.r_capbase=r_capbase;
    rates.r_delbase=r_delbase;
    rates.k_relbase=k_relbase;

    out_struct.rates=rates;
    out_struct.opts=Exp.opts;
    out_struct.fh1sizes=fh1sizes;
    out_struct.prmlocs=prmlocs;
    
end