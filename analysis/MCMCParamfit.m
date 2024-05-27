function params_all=MCMCParamfit(Exp,type)
% FORMINBAR  Makes and saves overview bargraphs for a collection of
% formins.
    %   out = MCMCParamfit(Exp) 
    %   
    %   Inputs:
    %         : 
    %
    %   
    %   
    %   See 

    arguments
        Exp
        type
    end

    NTMAX = 1e8;
    NTCHECK = 1000;
    KSCRITICAL = 0.01;
    NTADAPT = 100;
    
    [data,rates.k_capbase,rates.k_delbase,rates.r_capbase,rates.r_delbase,rates.k_relbase]=readinExp(Exp);
    
    % Pick initial guess 
    rng("shuffle")
    params=randi([1 10000],1,7);
    params(7)=rand;
    dx=ones(7);

    accepts=zeros(7);
    proposals=zeros(7);

    params_all=zeros(NTMAX,7);

    nt=0;
    currentntcheck=NTCHECK;
    [kpolys_nt,logll_nt]=loglikelihood(type,data,rates,proposals(1:6),proposals(7));
    while(nt<NTMAX)
        nt=nt+1;
        
       

        %randomly perturb one parameter
        [proposals,index]=generateproposal(params,dx);
        proposals(index)=proposals(index)+1;

        %calculate logll of proposal
        [kpolys_prop,logll_prop]=loglikelihood(type,data,rates,proposals(1:6),proposals(7));
        
        %Accept or reject proposal
        if logll_prop>logll_nt
            % Accept if new likelihood is greater
            params=proposals;
            logll_nt=logll_prop;
            kpolys_nt=kpolys_prop;
            accepts(index)=accepts(index)+1;
        elseif rand < exp(logll_prop-logll_nt)
            % Boltzmann test
            params=proposals;
            logll_nt=logll_prop;
            kpolys_nt=kpolys_prop;
            accepts(index)=accepts(index)+1;
        end
        params_all(nt,:)=params;

        if nt==3*currentntcheck
            ks=arrayfun(@(x) kstest2(x(currentntcheck:2*currentntcheck),x(2*currentntcheck:3*currentntcheck),'Alpha',KSCRITICAL),params_all); % will this work on matrix???
            if ks
                return
            else
                currentntcheck=currentntcheck*9;
            end
        end
    end
    
    % Still to do: at every NTAdpat before NTCHECK, see if p(accept)= 0.44, set
    % dx=dx(p(accept)/0.44)

end
function [data,k_capbase,k_delbase,r_capbase,r_delbase,k_relbase]=readinExp(Exp)
    % do something to read in formin info and stuff to take it out of
    % struct format
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
        % then you can just take the final values and times them by the
        % constants
        % and for rcap=C1*e^(-C2n), you can do rcap(C1,C2)=C1*(rcap(1,1))^C2

    % Return to original opts values 
    Exp.opts.k_cap=ogtable.k_cap;
    Exp.opts.k_del=ogtable.k_del;
    Exp.opts.r_cap=ogtable.r_cap;
    Exp.opts.r_del=ogtable.r_del;
    Exp.opts.k_rel=ogtable.k_rel;
end
function [kpolys,logll]=loglikelihood(type,data,rates,params,sigma)
    kpolys=calckpolys(type,rates,params);    
    
    data=struct2table(data);
    expdata=data.value;

    simdata=zeros(size(expdata));
    for i=1:length(expdata)
        if data.type(i)=="double"
            simdata(i)=kpolys.double(i);
        elseif data.type(i)=="ratio"
            simdata(i)=kpolys.dimer(i)/kpolys.double(i);
        elseif data.type(i)=="dimer"
            simdata(i)=kpolys.dimer(i);
        elseif data.type(i)=="single"
            simdata(i)=kpolys.single(i);
        else
            error('Error. \nNo valid experimental data type.')
        end
    end
    singleloglikelihoods=log(normpdf(expdata,simdata,sigma));
    logll=sum(singleloglikelihoods);
end

function out=calckpolys(type,rates,params)
    % calulcates log likelihood values for input parameters

    % calculate per PRM rates
    kcaps=cellfun(@(x) arrayfun(@(z) z*params(1),x), rates.k_capbase,'UniformOutput',false);
    kdels=cellfun(@(x) arrayfun(@(z) z*params(2),x), rates.k_delbase,'UniformOutput',false);
    rcaps=cellfun(@(x) arrayfun(@(z) ((z)^params(4))*params(3),x), rates.r_capbase,'UniformOutput',false);
    rdels=cellfun(@(x) arrayfun(@(z) z*params(5),x), rates.r_delbase,'UniformOutput',false);
    krels=cellfun(@(x) arrayfun(@(z) z*params(6),x), rates.k_relbase,'UniformOutput',false);
    kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) arrayfun(@(kcap,kdel,rcap,rdel,krel)kpolymerization(type,kcap,kdel,rcap,rdel,krel),kcap,kdel,rcap,rdel,krel),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins

    % sum up PRMs
    for i=1:length(kpolys)
        if length(kpolys{i})>1
            curkpoly=kpolys{i};
            kpolys{i}=curkpoly(1);
            for j=2:length(curkpoly)
                kpolys{i}=kpolys{i}+curkpoly(j);
            end
        end
    end

    % Sum up filaments
    kpolys=cellfun(@(obj) obj.add_fils,kpolys);

    out.single=arrayfun(@(x) x.single,kpolys);
    out.double=arrayfun(@(x) x.double,kpolys);
    out.dimer=arrayfun(@(x) x.dimer,kpolys);
end
function [proposals,i]=generateproposal(params,dx)
    i=randi([1 length(params)]);
    proposals=params;
    proposals(i)=proposals(i)+(dx(i)*(2*rand-1));
end