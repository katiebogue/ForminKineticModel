function params_out=MCMCParamfit(Exp,exptype,type,errtype,logtf,NTCHECK,NTADAPT,NTMAX,KSCRITICAL)
% 
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
        exptype % 1= is an experiment, 2= is a struct with rates, data, and resultsfolder and resultsdir
        type
        errtype % 1= use separate sigma values, 2= divide each by SEM
        logtf % whether or not to have step sizes in log space
        NTCHECK = 1000
        NTADAPT =50
        NTMAX =1e8
        KSCRITICAL =0.01
    end

    disp("MCMCParamfit.m starting") 

    % set up place to store data
    if exptype==1
        opts=Exp.opts;
        out_struct=readinExp(Exp);
        rates=out_struct.rates;
        data=out_struct.data;
        opts.update_results_folder
        opts.resultsfolder=strcat("MCMC_",opts.resultsfolder);
        clear Exp
        disp("Read in information from Experiment object")
    elseif exptype==2
        opts.resultsdir=Exp.resultsdir;
        opts.resultsfolder=strcat("MCMC_",Exp.resultsfolder);
        rates=Exp.rates;
        data=Exp.data;
        clear Exp
        disp("Loaded in pre-determined rates, data, and opts")
    else
        error("invalid exptype")
    end

    mkdir (opts.resultsdir,opts.resultsfolder)
    disp("created output directory")
    wkspc=fullfile(opts.resultsdir,opts.resultsfolder,"mcmc_results.mat");
    
    if errtype==1
        nsigma=length(unique(struct2table(data).type));
    elseif errtype==2
        nsigma=1;
    else
        error("invalid error type")
    end
    nparams=6+nsigma;

    % Pick initial guess 
    rng("shuffle")
    params=randi([1 10000],1,nparams);
    params(7:nparams)=rand(1,nparams-6);
    params(4)=randi([1 10],1,1);
    if logtf
        dx=zeros(nparams,1);
        dx(1:6)=[5,5,15,0,10,5];
    else
        dx=ones(nparams,1);
        dx(1:6)=[10^5,10^5,10^15,1,10^10,10^5];
    end
    
    accepts=zeros(nparams,1);
    proposals=zeros(nparams,1);

    nt=0;
    disp("saving workspace...")
    save(wkspc,'-v7.3')
    disp("saved workspace")
    m = matfile(wkspc,'Writable',true);
    disp("created matfile")
    m.params_all(NTMAX,nparams)=0;
    disp("created params_all matrix")
    params_temp=zeros(NTCHECK,nparams);

   
    nt_temp=0;
    last_nt=0;
    currentntcheck=NTCHECK;
    [kpolys_nt,logll_nt]=loglikelihood(type,data,rates,params(1:6),params(7:nparams),errtype);
    disp("starting MCMC loop")
    while(nt<NTMAX)
        nt=nt+1;
        nt_temp=nt_temp+1;
        
        if (nt<NTCHECK)&&(rem(nt,NTADAPT)==0)
            p_accept=accepts./proposals;
            p_accept(p_accept==0)=0.01;
            dx=dx.*p_accept/0.44;
            m.dx=dx;
        end

        %randomly perturb one parameter
        [proposal,index]=generateproposal(params,dx,logtf);
        proposals(index)=proposals(index)+1;

        %calculate logll of proposal
        [kpolys_prop,logll_prop]=loglikelihood(type,data,rates,proposal(1:6),proposal(7:nparams),errtype);
        
        %Accept or reject proposal
        if logll_prop>logll_nt
            % Accept if new likelihood is greater
            params=proposal;
            logll_nt=logll_prop;
            %kpolys_nt=kpolys_prop;
            accepts(index)=accepts(index)+1;
        elseif rand < exp(logll_prop-logll_nt)
            % Boltzmann test, Accept 
            params=proposal;
            logll_nt=logll_prop;
            %kpolys_nt=kpolys_prop;
            accepts(index)=accepts(index)+1;
        end
        params_temp(nt_temp,:)=params;

        if nt_temp==NTCHECK
            m.nt=nt;
            m.params_all(last_nt+1:nt,:)=params_temp;
            last_nt=nt;
            nt_temp=0;
            params_temp(:)=0;
            fprintf('appended params_all at nt: %d\n',nt)
        end

        if nt==3*currentntcheck
            ks=kstest2(m.params_all(currentntcheck:2*currentntcheck,3),m.params_all(2*currentntcheck:3*currentntcheck,3),'Alpha',KSCRITICAL); 
            if ~ks
                params_out=m.params_all(currentntcheck:3*currentntcheck,:);
                m.params_out=params_out;
                return
            else
                m.proposals=proposals;
                m.accepts=accepts;
                currentntcheck=currentntcheck*9;
                m.currentntcheck=currentntcheck;
                fprintf('nt: %d\n',nt)
            end
        end
    end
    params_out=0;
    save(m,'params_out')
end
function [kpolys,logll]=loglikelihood(type,data,rates,params,sigma,errtype)
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
    logll=0;
    if errtype==1
        types=unique(data.type);
        if ~length(types)==length(sigma)
            error("Number of unique data types not the same as number of sigmas")
        end
        for i=1:length(sigma)
            rows = matches(data.type,types(i));
            SSE=sum((expdata(rows)-simdata(rows)).^2);
            sigma2=sigma(i)^2;
            logll_i=(-SSE/(2*sigma2))-(0.5*length(expdata(rows))*log(2*pi*sigma2));
            logll=logll+logll_i;
        end
    elseif errtype==2
        for i=1:length(expdata)
            SEM=(data.errtop(i)+data.errbot(i))/2;
            SSE=((expdata(i)-simdata(i)).^2)/SEM;
            sigma2=sigma^2;
            logll_i=(-SSE/(2*sigma2))-(0.5*log(2*pi*sigma2));
            logll=logll+logll_i;
        end
    else
        error("invalid error type")
    end
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
function [proposals,i]=generateproposal(params,dx,logtf)
    i=randi([1 length(params)]);
    proposals=params;
    if logtf
        proposals(i)=proposals(i)+(10^(dx(i))*(2*rand-1));
    else
        proposals(i)=proposals(i)+(dx(i)*(2*rand-1));
    end
end