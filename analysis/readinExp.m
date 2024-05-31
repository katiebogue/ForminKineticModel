function out_struct=readinExp(Exp)
    % read in formin info 
    ogtable=struct('k_cap',Exp.opts.k_cap,'k_del',Exp.opts.k_del,'r_cap',Exp.opts.r_cap,'r_del',Exp.opts.r_del,'k_rel',Exp.opts.k_rel);
    Exp.opts.k_cap= 1;
    Exp.opts.k_del=1; 
    Exp.opts.r_cap=1; 
    Exp.opts.r_del=1; 
    Exp.opts.k_rel=1;

    data=Exp.data;
    L=length(data);

    k_capbase=cell([L 1]);
    k_delbase=cell([L 1]);
    r_capbase=cell([L 1]);
    r_delbase=cell([L 1]);
    k_relbase=cell([L 1]);

    for i=1:length(data)
        formin=data(i).formin;
        kcaps=FilType.empty(formin.PRMCount,0);
        kdels=FilType.empty(formin.PRMCount,0);
        rcaps=FilType.empty(formin.PRMCount,0);
        rdels=FilType.empty(formin.PRMCount,0);
        krels=FilType.empty(formin.PRMCount,0);
        for j=1:formin.PRMCount
            prm1=formin.PRMList(1,j);
            kcaps(j)=prm1.kcap;
            kdels(j)=prm1.kdel;
            rcaps(j)=prm1.rcap;
            rdels(j)=prm1.rdel;
            krels(j)=prm1.krel;
        end
        k_capbase{i,1}=kcaps;
        k_delbase{i,1}=kdels;
        r_capbase{i,1}=rcaps;
        r_delbase{i,1}=rdels;
        k_relbase{i,1}=krels;
    end
    % get values with all parameters set to 1

    % Return to original opts values 
    Exp.opts.k_cap=ogtable.k_cap;
    Exp.opts.k_del=ogtable.k_del;
    Exp.opts.r_cap=ogtable.r_cap;
    Exp.opts.r_del=ogtable.r_del;
    Exp.opts.k_rel=ogtable.k_rel;

    out_struct.data=data;

    rates.k_capbase=k_capbase;
    rates.k_delbase=k_delbase;
    rates.r_capbase=r_capbase;
    rates.r_delbase=r_delbase;
    rates.k_relbase=k_relbase;

    out_struct.rates=rates;
    out_struct.opts=Exp.opts;
end