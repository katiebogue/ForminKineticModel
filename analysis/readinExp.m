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
        kcaps=zeros(formin.PRMCount,5);
        kdels=zeros(formin.PRMCount,5);
        rcaps=zeros(formin.PRMCount,5);
        rdels=zeros(formin.PRMCount,5);
        krels=zeros(formin.PRMCount,5);
        for j=1:formin.PRMCount
            prm1=formin.PRMList(1,j);
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
    end
    % get values with all parameters set to 1

    % Return to original opts values 
    Exp.opts.k_cap=ogtable.k_cap;
    Exp.opts.k_del=ogtable.k_del;
    Exp.opts.r_cap=ogtable.r_cap;
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
end