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
        logtf % whether or not to have params in log space
        NTCHECK = 1000
        NTADAPT =100
        NTMAX =10^6
        KSCRITICAL =0.01
    end

    NBINS = 200;
    PARAMMAX = 20; % in log-space
    PARAMMIN = -5; % in log-space

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
    params(2)=0.1;
    params(4)=1;
    %params(4)=randi([1 5],1,1);
    if logtf
        params=log10(params);
        dx=ones(nparams,1)./100;
        dx(1:3)=[5,5,5];
        if type=="4st"
            dx(5:6)=[10,5];
        end
    else
        PARAMMIN=10^PARAMMIN;
        PARAMMAX=10^PARAMMAX;
        dx=ones(nparams,1);
        dx(1:3)=[10^5,10^5,10^15];
        if type=="4st"
            dx(5:6)=[10^10,10^5];
        end
    end
    
    accepts=zeros(nparams,1);
    proposals=zeros(nparams,1);
    accepts_temp=zeros(nparams,1);
    proposals_temp=zeros(nparams,1);
    accepts_adapt=zeros(nparams,1);
    proposals_adapt=zeros(nparams,1);

    nt=0;
    disp("saving workspace...")
    save(wkspc,'-v7.3')
    disp("saved workspace")
    m = matfile(wkspc,'Writable',true);
    disp("created matfile")
    m.params_all(NTMAX,nparams)=0;
    m.logparams_all(NTMAX,nparams)=0;
    m.paccept_matrix(ceil(NTMAX/NTCHECK)+1,nparams+1)=0;
    disp("created paccept_matrix")
    m.paramHistCounts_matrices(nparams,NBINS,ceil(NTMAX/NTCHECK))=0;
    disp("created paramHistCounts_matrices")
    m.paramHistCounts_matrices_nts(ceil(NTMAX/NTCHECK),1)=0;
    m.ksvals(ceil(NTMAX/NTCHECK),nparams)=0;
    disp("created matfile matrices")
    params_temp=zeros(NTCHECK,nparams);

   
    nt_temp=0;
    last_nt=0;
    ntcheck_count=0;
    ksvalindex=1;
    currentntcheck=NTCHECK;
    paramHistCounts = zeros(nparams,NBINS);
    paramHistCountsPrevious = zeros(nparams,NBINS);
    logll_nt=loglikelihood(type,data,rates,params(1:nkpolyparams),params(nkpolyparams+1:nparams),errtype,logtf);
    minlogll=logll_nt;
    disp("starting MCMC loop")
    while(nt<NTMAX)
        nt=nt+1;
        nt_temp=nt_temp+1;
        
        if (nt<NTCHECK-NTADAPT)&&(rem(nt,NTADAPT)==0)

            disp("adapting step size...")
            fprintf('nt: %d\n',nt)
            disp("total proposal count:")
            disp(proposals)
            disp("total acceptance count:")
            disp(accepts)

            proposals_new=proposals-proposals_adapt;
            accepts_new=accepts-accepts_adapt;

            disp("new proposal count:")
            disp(proposals_new)
            disp("new acceptance count:")
            disp(accepts_new)
            
            p_accept=accepts./proposals;
            disp("total acceptance probability:")
            disp(p_accept)
            p_accept=accepts_new./proposals_new;
            disp("new acceptance probability:")
            disp(p_accept)
            p_accept(p_accept==0)=0.01;
            p_accept(proposals==0)=0.44;
            in=(p_accept>0.6 | p_accept<0.3);
            disp("previous step sizes: ")
            disp(dx)
            dx(in)=dx(in).*(p_accept(in)./0.44);

            disp("new step sizes: ")
            disp(dx)
            m.dx=dx;
            proposals_adapt=proposals;
            accepts_adapt=accepts;
        end

        %randomly perturb one parameter
        [proposal,index]=generateproposal(params,dx,logtf);
        proposals(index)=proposals(index)+1;
        proposals_temp(index)=proposals_temp(index)+1;

        if proposal(index)<PARAMMIN || proposal(index)>PARAMMAX
            % reject anything beyond the boundaries
            params_temp(nt_temp,:)=params;
        else
            %calculate logll of proposal
            logll_prop=loglikelihood(type,data,rates,proposal(1:nkpolyparams),proposal(nkpolyparams+1:nparams),errtype,logtf);
            
            %Accept or reject proposal
            if logll_prop>logll_nt
                % Accept if new likelihood is greater
                params=proposal;
                logll_nt=logll_prop;
                %kpolys_nt=kpolys_prop;
                accepts(index)=accepts(index)+1;
                accepts_temp(index)=accepts_temp(index)+1;
                if logll_nt>minlogll
                    minlogll=logll_nt;
                    minlogll_params=params;
                end
            elseif rand < exp(logll_prop-logll_nt)
                % Boltzmann test, Accept 
                params=proposal;
                logll_nt=logll_prop;
                %kpolys_nt=kpolys_prop;
                accepts(index)=accepts(index)+1;
                accepts_temp(index)=accepts_temp(index)+1;
            end
            params_temp(nt_temp,:)=params;
        end

         if nt==NTCHECK
            binEdges=zeros(nparams,NBINS+1);
            Y = zeros(nparams,NBINS);
            for i=1:nparams
                %maxp=max(params_temp(:,i));
                %minp=min(params_temp(:,i));
                [Y(i,:),binEdges(i,:)]=histcounts(params_temp((NTCHECK/2):end,i),NBINS);
            end
            m.paramHistCounts_matrices(:,:,1)=Y;
            m.paramHistCounts_matrices_nts(1,1)=nt;
            m.binEdges=binEdges;
            clear Y
            HistCountIndex=2;
        end

        if nt_temp==NTCHECK
            ntcheck_count=ntcheck_count+1;
            m.nt=nt;
            if logtf
                m.logparams_all(last_nt+1:nt,:)=params_temp;
                params_temp=10.^params_temp;
            end
            m.params_all(last_nt+1:nt,:)=params_temp;
            m.minlogll=minlogll;
            m.minlogll_params=minlogll_params;
            m.paccept_matrix(ntcheck_count,1)=nt;
            m.paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            accepts_temp(:)=0;
            proposals_temp(:)=0;
            last_nt=nt;
            nt_temp=0;
            params_temp(:)=0;
            fprintf('appended params_all at nt: %d\n',nt)
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
            fprintf('nt: %d\n',nt)
            fprintf('currentntcheck: %d\n',currentntcheck)
            m.paramHistCounts_matrices_nts(HistCountIndex,1)=currentntcheck*2;
            m.paramHistCounts_matrices_nts(HistCountIndex+1,1)=currentntcheck*3;
            m.paramHistCounts_matrices(:,:,HistCountIndex)=paramHistCountsPrevious;
            m.paramHistCounts_matrices(:,:,HistCountIndex+1)=paramHistCounts;
            HistCountIndex=HistCountIndex+2;
            
            disp('Begining KS test')
            % intialize
            ksStatistic = zeros(nparams,1);
            cdf1 = zeros(nparams,NBINS);
            cdf2 = zeros(nparams,NBINS);
            for i=1:nparams
                % compute cumulative distribution functions
                cdf1(i,:) = cumsum(paramHistCountsPrevious(i,:))./(currentntcheck);
                cdf2(i,:) = cumsum(paramHistCounts(i,:))./(currentntcheck);
                % compute ksStatistic
                ksStatistic(i) = max(abs(cdf1(i,:)-cdf2(i,:)));
                %[ks(i),p,ks2stat]=kstest2(paramHistCountsPrevious(i,:)./currentntcheck,paramHistCounts(i,:)./currentntcheck,'Alpha',KSCRITICAL); 
                %fprintf('%d ks: %d p: %d ks2stat: %d\n',i,ks(i),p,ks2stat)
            end
            m.ksvals(ksvalindex,:)=ksStatistic';
            ksvalindex=ksvalindex+1;
            %ks=kstest2(m.params_all(currentntcheck:2*currentntcheck,3),m.params_all(2*currentntcheck:3*currentntcheck,3),'Alpha',KSCRITICAL); 
            if all(ksStatistic(:) < KSCRITICAL)
                disp('KS test successful')
                m.proposals=proposals;
                m.accepts=accepts;
                m.ksStatistic=ksStatistic;
                m.params_out=m.params_all(currentntcheck:3*currentntcheck,:);
                m.paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
                m.paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
                %visualizePosteriors(fullfile(opts.resultsdir,opts.resultsfolder),1)
                return
            else
                disp(ksStatistic)
                disp('KS test unsuccessful')
                m.proposals=proposals;
                m.accepts=accepts;
                currentntcheck=currentntcheck*3;
                m.currentntcheck=currentntcheck;
                fprintf('new currentntcheck: %d\n',currentntcheck)
                paramHistCounts = zeros(nparams,NBINS);
                paramHistCountsPrevious = zeros(nparams,NBINS);
                fprintf('nt: %d\n',nt)
            end
        end
    end
    m.paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
    m.paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
    m.params_out=0;
    visualizePosteriors(fullfile(opts.resultsdir,opts.resultsfolder),1)
end
function logll=loglikelihood(type,data,rates,params,sigma,errtype,logtf)
    SSE=sum((1.5-params).^2);
    sigma2=sigma(1)^2;
    logll=(-SSE/(2*sigma2))-(0.5*length(params)*log(2*pi*sigma2));
    return

    if logtf
        kpolys=calckpolys(type,rates,10.^(params)); 
    else
        kpolys=calckpolys(type,rates,params); 
    end
    
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
    proposals(i)=proposals(i)+(dx(i)*(2*rand-1));
    % if logtf
    %     proposals(i)=proposals(i)+(10^(dx(i))*(2*rand-1));
    % else
    %     proposals(i)=proposals(i)+(dx(i)*(2*rand-1));
    % end
end