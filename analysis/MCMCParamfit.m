function MCMCParamfit(Exp,exptype,type,errtype,logtf,NTCHECK,NTADAPT,NTMAX,KSCRITICAL)
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
        NTMAX =10000
        KSCRITICAL =0.01
    end

    NBINS = 500;

    disp("MCMCParamfit.m starting") 

    % set up place to store data
    if exptype==1
        opts=Exp.opts;
        out_struct=readinExp(Exp);
        rates=out_struct.rates;
        datatab=out_struct.data;
        opts.update_results_folder
        opts.resultsfolder=strcat("MCMC_",opts.resultsfolder);
        clear Exp
        disp("Read in information from Experiment object")
    elseif exptype==2
        opts.resultsdir=Exp.resultsdir;
        opts.resultsfolder=strcat("MCMC_",Exp.resultsfolder);
        rates=Exp.rates;
        datatab=Exp.data;
        clear Exp
        disp("Loaded in pre-determined rates, data, and opts")
    else
        error("invalid exptype")
    end

    data=struct2table(datatab);

    mkdir (opts.resultsdir,opts.resultsfolder)
    disp("created output directory")
    wkspc=fullfile(opts.resultsdir,opts.resultsfolder,"mcmc_results.mat");
    
    if errtype==1
        nsigma=length(unique(data.type));
    elseif errtype==2
        nsigma=1;
    else
        error("invalid error type")
    end
    if type=="3st"
        nkpolyparams=4;
    elseif type=="4st"
        nkpolyparams=6;
    end
    nparams=nkpolyparams+nsigma;

    % Pick initial guess 
    rng("shuffle")
    params=randi([1 1000],1,nparams);
    params(nkpolyparams+1:nparams)=rand(1,nparams-nkpolyparams);
    params(4)=1;
    %params(4)=randi([1 5],1,1);
    if logtf
        dx=ones(nparams,1)./100;
        dx(1:3)=[5,5,5];
        if type=="4st"
            dx(5:6)=[10,5];
        end
    else
        dx=ones(nparams,1);
        dx(1:3)=[10^5,10^5,10^15];
        if type=="4st"
            dx(5:6)=[10^10,10^5];
        end
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
    paramHistCounts = zeros(nparams,NBINS);
    paramHistCountsPrevious = zeros(nparams,NBINS);
    logll_nt=loglikelihood(type,data,rates,params(1:nkpolyparams),params(nkpolyparams+1:nparams),errtype);
    disp("starting MCMC loop")
    while(nt<NTMAX)
        nt=nt+1;
        nt_temp=nt_temp+1;
        
        if (nt<NTCHECK)&&(rem(nt,NTADAPT)==0)
            p_accept=accepts./proposals;
            p_accept(p_accept==0)=0.01;
            p_accept(proposals==0)=0.44;
            in=(dx>0.58 | dx<0.3);
            dx(in)=dx(in).*(p_accept(in)./0.44);
            m.dx=dx;
        end

        %randomly perturb one parameter
        [proposal,index]=generateproposal(params,dx,logtf);
        proposals(index)=proposals(index)+1;

        %dont allow negative values
        proposal(proposal<0)=eps(0);

        %calculate logll of proposal
        logll_prop=loglikelihood(type,data,rates,proposal(1:nkpolyparams),proposal(nkpolyparams+1:nparams),errtype);
        
        %Accept or reject proposal
        if logll_prop>logll_nt
            % Accept if new likelihood is greater
            params=proposal;
            logll_nt=logll_prop;
            %kpolys_nt=kpolys_prop;
            accepts(index)=accepts(index)+1;
            minlogll=logll_nt;
            minlogll_params=params;
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
            m.minlogll=minlogll;
            m.minlogll_params=minlogll_params;
            last_nt=nt;
            nt_temp=0;
            params_temp(:)=0;
            fprintf('appended params_all at nt: %d\n',nt)
        end

        if nt==NTCHECK
            binEdges=zeros(nparams,NBINS+1);
            for i=1:nparams
                %maxp=max(params_temp(:,i));
                %minp=min(params_temp(:,i));
                [Y,binEdges(i,:)]=discretize(params_temp(:,i),NBINS);
            end
            clear Y
        end
        if nt<=currentntcheck
        elseif nt<=2*currentntcheck
            for i=1:nparams
                if params(i)<binEdges(i,1)
                    index=1;
                elseif params(i)>binEdges(i,NBINS+1)
                    index=NBINS;
                else
                    index=discretize(params(i),binEdges(i,:));
                end
                paramHistCountsPrevious(i,index)=paramHistCountsPrevious(i,index)+1;
            end
            %disp('Updated bins for 2nd segment')
        elseif nt<=3*currentntcheck
            for i=1:nparams
                if params(i)<binEdges(i,1)
                    index=1;
                elseif params(i)>binEdges(i,NBINS+1)
                    index=NBINS;
                else
                    index=discretize(params(i),binEdges(i,:));
                end
                paramHistCounts(i,index)=paramHistCounts(i,index)+1;
            end
            %disp('Updated bins for 3rd segment')
        end

        if nt==3*currentntcheck
            disp('Begining KS test')
            ks=zeros(nparams,1);
            for i=1:nparams
                [ks(i),p,ks2stat]=kstest2(paramHistCountsPrevious(i,:)./currentntcheck,paramHistCounts(i,:)./currentntcheck,'Alpha',KSCRITICAL); 
                fprintf('%d ks: %d p: %d ks2stat: %d\n',i,ks(i),p,ks2stat)
            end
            %ks=kstest2(m.params_all(currentntcheck:2*currentntcheck,3),m.params_all(2*currentntcheck:3*currentntcheck,3),'Alpha',KSCRITICAL); 
            if ~all(ks)
                disp('KS test successful')
                m.proposals=proposals;
                m.accepts=accepts;
                m.params_out=m.params_all(currentntcheck:3*currentntcheck,:);
                %visualizePosteriors(opts.resultsdir,1)
                return
            else
                disp('KS test unsuccessful')
                m.proposals=proposals;
                m.accepts=accepts;
                currentntcheck=currentntcheck*9;
                m.currentntcheck=currentntcheck;
                paramHistCounts = zeros(n_params,NBINS);
                paramHistCountsPrevious = zeros(n_params,NBINS);
                fprintf('nt: %d\n',nt)
            end
        end
    end
    m.params_out=0;
    %visualizePosteriors(opts.resultsdir,1)
end
function logll=loglikelihood(type,data,rates,params,sigma,errtype)
    kpolys=calckpolys(type,rates,params);    
    
    expdata=data.value;

    simdata=zeros(size(expdata));
    for i=1:length(expdata)
        if data.type(i)=="double"
            simdata(i)=kpolys{i}(1,2);
        elseif data.type(i)=="ratio"
            simdata(i)=kpolys{i}(1,3)/kpolys{i}(1,2);
        elseif data.type(i)=="dimer"
            simdata(i)=kpolys{i}(1,3);
        elseif data.type(i)=="single"
            simdata(i)=kpolys{i}(1,1);
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

function kpolys=calckpolys(type,rates,params)
    % calulcates log likelihood values for input parameters

    % calculate per PRM rates
    kcaps=cellfun(@(x) x.*params(1), rates.k_capbase,'UniformOutput',false);
    kdels=cellfun(@(x) x.*params(2), rates.k_delbase,'UniformOutput',false);
    rcaps=cellfun(@(x) ((x).^params(4)).*params(3), rates.r_capbase,'UniformOutput',false);
    if type=="4st"
        rdels=cellfun(@(x) x.*params(5), rates.r_delbase,'UniformOutput',false);
        krels=cellfun(@(x) x.*params(6), rates.k_relbase,'UniformOutput',false);
    elseif type=="3st"
        rdels=kcaps;
        krels=kcaps;
    end
    kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) kpolymerization(type,kcap,kdel,rcap,rdel,krel),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins

    for i=1:length(kpolys)
        PRMsum=sum(kpolys{i},1); % sum up PRMs
        kpolys{i}=[PRMsum(1),sum(PRMsum(2:3)),sum(PRMsum(4:5))]; % Sum up filaments
    end
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