function MCMCParamfit(Exp,exptype,type,errtype,matfileTF, NTCHECK,NTADAPT,NTMAX,KSCRITICAL,nondim,prcalc)
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
        exptype % 1= is an experiment, 2= is a struct with rates, data, and resultsfolder and resultsdir, and fh1sizes and prmlocs
        type
        errtype % 1= use separate sigma values, 2= divide each by SEM
        matfileTF=0 % whether to save to a matfile or to keep things in memory
        NTCHECK = 3000
        NTADAPT =100
        NTMAX =10^7
        KSCRITICAL =0.02
        nondim= 0 % whether to use nondimensionality
        prcalc= 1
    end

    
    NBINS = 200;
    PARAMMAX = 20; % in log-space
    PARAMMIN = -10; % in log-space
    SIGMAMAX = 2; % in log-space
    SIGMAMIN = -2; % in log-space
    EXPMIN = 0.1; % non log-space, applies to the 4th parameter (or 3rd if nondimensional)
    EXPMAX = 10; % non log-space, applies to the 4th parameter (or 3rd if nondimensional)
    XLOCMAX= 35.5;
    XLOCMIN= 0;
    YLOCMAX= 20;
    YLOCMIN= 0;

    disp("MCMCParamfit.m starting") 

    % set up place to store data
    if exptype==1
        opts=Exp.opts;
        if prcalc
            opts.set_equation(2);
        end
        out_struct=readinExp(Exp);
        opts.set_equation(1);
        rates=out_struct.rates;
        datatab=out_struct.data;
        opts.update_results_folder
        opts.resultsfolder=strcat("MCMC_",opts.resultsfolder);
        fh1lengths=out_struct.fh1sizes;
        prmlocs=out_struct.prmlocs;
        clear Exp
        disp("Read in information from Experiment object")
    elseif exptype==2
        opts.resultsdir=Exp.resultsdir;
        opts.resultsfolder=strcat("MCMC_",Exp.resultsfolder);
        rates=Exp.rates;
        datatab=Exp.data;
        fh1lengths=Exp.fh1sizes;
        prmlocs=Exp.prmlocs;
        clear Exp
        disp("Loaded in pre-determined rates, data, and opts")
    else
        error("invalid exptype")
    end

    data=struct2table(datatab);
    divdatapoint=0;
    if nondim
        for a=1:height(data)
            if length(data.groups{a})==2
                if ~divdatapoint
                    divdatapoint=a;
                else
                    error("two possible datapoints to divide by")
                end
            end
        end
        if ~divdatapoint
            error("no datapoint matches condition to be the dividing point")
        end
        divkpoly=data.value(divdatapoint);
        data.ratiovalues=data.value;
        for a=1:height(data)
            if data.type(a)=="double"
                data.ratiovalues(a)=data.value(a)/divkpoly;
            end
        end
    end

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
        if nondim
            nkpolyparams=3;
        else
            nkpolyparams=4;
        end
    elseif type=="4st"
        if nondim
            nkpolyparams=5;
        else
            nkpolyparams=6;
        end
    end
    if prcalc
        nkpolyparams=nkpolyparams+2;
    end
    nparams=nkpolyparams+nsigma;

    % Pick initial guess
    % note: all in log10 space
    rng("shuffle")
    % params=rand([1 nparams])+1;
    % dx=ones(nparams,1);
    params=zeros(1,nparams);
    params(nkpolyparams+1:nparams)=rand(1,nparams-nkpolyparams); % set initial sigma values to between 0-1
    dx=ones(nparams,1);
    dx(1:3)=[0.1,0.1,0.1];
    if nondim
        params(1)=randi([11 14],1,1); % alpha_del initial
        params(2)=-rand; % beta_cap initial
        params(3)=1; %rcap exp initial
        dx(3)=1;
    else
        params(1)=randi([11 14],1,1); % kcap initial
        params(2)=-rand; % kdel initial
        params(3)=randi([13 15],1,1); %rcap initial
        params(4)=1; %rcap exp initial
        dx(4)=1;
    end
    if type=="4st"
        if nondim
            params(4)=1; %gamma_del initial
            params(5)=randi([3 5],1,1); %tau_rel initial
            dx(5:6)=[1,5];
        else
            params(5)=randi([10 15],1,1); %rdel initial
            params(6)=randi([3 5],1,1); %krel initial
            dx(5:6)=[5,5];
        end
    end

    if prcalc
        params(nkpolyparams-1:nkpolyparams)=0; %initial delivery location
    end

    accepts=zeros(nparams,1);
    proposals=zeros(nparams,1);
    accepts_temp=zeros(nparams,1);
    proposals_temp=zeros(nparams,1);
    accepts_adapt=zeros(nparams,1);
    proposals_adapt=zeros(nparams,1);

    step_up=zeros(nparams,1);
    step_down=zeros(nparams,1);

    nt=0;

    if matfileTF
        disp("saving workspace...")
        save(wkspc,'-v7.3')
        disp("saved workspace")
        m = matfile(wkspc,'Writable',true);
        disp("created matfile")
        m.logparams_all(NTMAX,nparams)=0;
        m.paccept_matrix(ceil(NTMAX/NTCHECK)+1,nparams+1)=0;
        disp("created paccept_matrix")
        m.paramHistCounts_matrices(nparams,NBINS,ceil(NTMAX/NTCHECK))=0;
        disp("created paramHistCounts_matrices")
        m.paramHistCounts_matrices_nts(ceil(NTMAX/NTCHECK),1)=0;
        m.ksvals(ceil(NTMAX/NTCHECK),nparams)=0;
        disp("created matfile matrices")
        params_temp=zeros(3*NTCHECK,nparams);
    else
        logparams_all(NTMAX,nparams)=0;
        paccept_matrix(ceil(NTMAX/NTCHECK)+1,nparams+1)=0;
        paramHistCounts_matrices(nparams,NBINS,ceil(NTMAX/NTCHECK))=0;
        paramHistCounts_matrices_nts(ceil(NTMAX/NTCHECK),1)=0;
        ksvals(ceil(NTMAX/NTCHECK),nparams)=0;
    end

   
    nt_temp=0;
    last_nt=0;
    ntcheck_count=0;
    ksvalindex=1;
    currentntcheck=NTCHECK;
    paramHistCounts = zeros(nparams,NBINS);
    paramHistCountsPrevious = zeros(nparams,NBINS);
    [logll_nt,divvalue]=loglikelihood(type,data,rates,params(1:nkpolyparams),params(nkpolyparams+1:nparams),errtype,nondim,divdatapoint,prcalc,prmlocs,fh1lengths);
    minlogll=logll_nt;
    
    disp("starting MCMC loop")
    while(nt<NTMAX)
        nt=nt+1;
        nt_temp=nt_temp+1;
        
        if (nt<NTCHECK-NTADAPT)&&(rem(nt,NTADAPT)==0)

            disp("adapting step size...")
            fprintf('nt: %d\n',nt)
            % disp("total proposal count:")
            % disp(proposals')
            % disp("total acceptance count:")
            % disp(accepts')

            proposals_new=proposals-proposals_adapt;
            accepts_new=accepts-accepts_adapt;

            % disp("new proposal count:")
            % disp(proposals_new')
            % disp("new acceptance count:")
            % disp(accepts_new')
            
            p_accept=accepts./proposals;
            % disp("total acceptance probability:")
            % disp(p_accept')
            p_accept=accepts_new./proposals_new;
            disp("new acceptance probability:")
            disp(p_accept')
            p_accept(p_accept==0)=0.01;
            p_accept(proposals==0)=0.44;
            in=(p_accept>0.6 | p_accept<0.3);
            disp("previous step sizes: ")
            disp(dx')
            dx(in)=dx(in).*(p_accept(in)./0.44);

            disp("new step sizes: ")
            disp(dx')
            if matfileTF
                m.dx=dx;
            end
            proposals_adapt=proposals;
            accepts_adapt=accepts;
        end

        %randomly perturb one parameter
        [proposal,index]=generateproposal(params,dx);
        proposals(index)=proposals(index)+1;
        proposals_temp(index)=proposals_temp(index)+1;

        if (prcalc && index==(nkpolyparams) && (proposal(index)>YLOCMAX || proposal(index)<YLOCMIN))... % del location thresholds
                || (prcalc && index==(nkpolyparams-1) && (proposal(index)>XLOCMAX || proposal(index)<XLOCMIN))... %del location thresholds
                || (index~=(4-nondim) && proposal(index)<PARAMMIN && index<=(nkpolyparams-prcalc-prcalc)) ...% regular param thresholds
                || (index~=(4-nondim) && proposal(index)>PARAMMAX && index<=(nkpolyparams-prcalc-prcalc)) ... % regular param thresholds
                || (index>nkpolyparams && (proposal(index)>SIGMAMAX || proposal(index)<SIGMAMIN)) ... % sigma thresholds
                || (index==(4-nondim) && (proposal(index)>EXPMAX || proposal(index)<EXPMIN)) % rcap exponent threshold
            
            % reject anything beyond the boundaries
            params_temp(nt_temp,:)=params;
        else
            %calculate logll of proposal
            [logll_prop,divvalue]=loglikelihood(type,data,rates,proposal(1:nkpolyparams),proposal(nkpolyparams+1:nparams),errtype,nondim,divdatapoint,prcalc,prmlocs,fh1lengths);
            % trueparams=gettrueparams(proposal,divvalue,divkpoly,nkpolyparams);
            % exptrueparams=10.^(trueparams);
            % exptrueparams(4)=trueparams(4);
            % kps=calckpolys(type,rates,exptrueparams,0);
            % val=kps{divdatapoint}(1,2);
            % if abs(val-divkpoly)>10^-5
            %     error("not getting to correct kpoly")
            % end
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
                    if nondim
                        minlogll_params=gettrueparams(params,divvalue,divkpoly,nkpolyparams,prcalc);
                        minlogll_params_raw=params;
                    else
                        minlogll_params=params;
                    end
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
            if matfileTF
                m.paramHistCounts_matrices(:,:,1)=Y;
                m.paramHistCounts_matrices_nts(1,1)=nt;
                m.binEdges=binEdges;
            else
                paramHistCounts_matrices(:,:,1)=Y;
                paramHistCounts_matrices_nts(1,1)=nt;
            end
            clear Y
            HistCountIndex=2;
        end

        if nt_temp==3*NTCHECK
            ntcheck_count=ntcheck_count+1;
            if matfileTF
                m.nt=nt;
                m.logparams_all(last_nt+1:nt,:)=params_temp;
                m.minlogll=minlogll;
                m.minlogll_params=minlogll_params;
                if nondim
                    m.minlogll_params_raw=minlogll_params_raw;
                end
                m.paccept_matrix(ntcheck_count,1)=nt;
                m.paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            else
                logparams_all(last_nt+1:nt,:)=params_temp;
                paccept_matrix(ntcheck_count,1)=nt;
                paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            end
            disp([accepts_temp./proposals_temp]')
            accepts_temp(:)=0;
            proposals_temp(:)=0;
            last_nt=nt;
            nt_temp=0;
            params_temp(:)=0;
            % disp("step up:")
            % disp(step_up')
            % disp("step down:")
            % disp(step_down')
            fprintf('appended logparams_all at nt: %d\n',nt)
            %makeupdateplot
        end

        if nt==NTMAX
            ntcheck_count=ntcheck_count+1;
            if matfileTF
                m.nt=nt;
                m.logparams_all(last_nt+1:nt,:)=params_temp(1:nt_temp,:);
                m.minlogll=minlogll;
                m.minlogll_params=minlogll_params;
                if nondim
                    m.minlogll_params_raw=minlogll_params_raw;
                end
                m.paccept_matrix(ntcheck_count,1)=nt;
                m.paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            else
                logparams_all(last_nt+1:nt,:)=params_temp(1:nt_temp,:);;
                paccept_matrix(ntcheck_count,1)=nt;
                paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            end
            disp([accepts_temp./proposals_temp]')
            accepts_temp(:)=0;
            proposals_temp(:)=0;
            last_nt=nt;
            nt_temp=0;
            params_temp(:)=0;
            % disp("step up:")
            % disp(step_up')
            % disp("step down:")
            % disp(step_down')
            fprintf('appended logparams_all at nt: %d\n',nt)
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
            if matfileTF
                m.paramHistCounts_matrices_nts(HistCountIndex,1)=currentntcheck*2;
                m.paramHistCounts_matrices_nts(HistCountIndex+1,1)=currentntcheck*3;
                m.paramHistCounts_matrices(:,:,HistCountIndex)=paramHistCountsPrevious;
                m.paramHistCounts_matrices(:,:,HistCountIndex+1)=paramHistCounts;
            else
                paramHistCounts_matrices_nts(HistCountIndex,1)=currentntcheck*2;
                paramHistCounts_matrices_nts(HistCountIndex+1,1)=currentntcheck*3;
                paramHistCounts_matrices(:,:,HistCountIndex)=paramHistCountsPrevious;
                paramHistCounts_matrices(:,:,HistCountIndex+1)=paramHistCounts;
            end
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
                
            end
            if matfileTF
                m.ksvals(ksvalindex,:)=ksStatistic';
            else
                ksvals(ksvalindex,:)=ksStatistic';
            end
            ksvalindex=ksvalindex+1;
            if all(ksStatistic(:) < KSCRITICAL)
                disp('KS test successful')
                if matfileTF
                    m.proposals=proposals;
                    m.accepts=accepts;
                    m.ksStatistic=ksStatistic;
                    m.parameters_all=m.logparams_all(currentntcheck:3*currentntcheck,:);
                    m.paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
                    m.paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
                else
                    parameters_all=logparams_all(currentntcheck:3*currentntcheck,:);
                    paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
                    paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
                    save(wkspc)
                    disp("saved workspace to .mat")
                end
                visualizePosteriors(fullfile(opts.resultsdir,opts.resultsfolder),1)
                return
            else
                disp(ksStatistic)
                disp('KS test unsuccessful')
                currentntcheck=currentntcheck*3;
                if matfileTF
                    m.currentntcheck=currentntcheck;
                    m.proposals=proposals;
                    m.accepts=accepts;
                end
                fprintf('new currentntcheck: %d\n',currentntcheck)
                paramHistCounts = zeros(nparams,NBINS);
                paramHistCountsPrevious = zeros(nparams,NBINS);
                %fprintf('nt: %d\n',nt)
                %makeupdateplot
                save(wkspc)
                disp("saved workspace to .mat")
            end
        end
    end
    if matfileTF
        m.paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
        m.paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
        m.parameters_all=0;
    else
        paramHistCounts_matrices(:,:,HistCountIndex:end)=[];
        paramHistCounts_matrices_nts(HistCountIndex:end,:)=[];
        parameters_all=0;
        save(wkspc)
        disp("saved workspace to .mat")
    end
    visualizePosteriors(fullfile(opts.resultsdir,opts.resultsfolder),1)

    function [proposals,i]=generateproposal(params,dx)
        i=randi([1 length(params)]);
        proposals=params;
        step=dx(i)*(2*rand-1);
        if step>0
            step_up(i)=step_up(i)+1;
        else
            step_down(i)=step_down(i)+1;
        end
        proposals(i)=proposals(i)+step;
    end

    function makeupdateplot
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs iteration');hold on;
        for k=1:nparams
            subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p);
            if matfileTF
                plot(m.logparams_all(1:nt,k),'LineWidth',0.05);
            else
                plot(logparams_all(1:nt,k),'LineWidth',0.05);
            end
            p=p+1;
        end
        pause(20)
    end

end


function [logll,divvalue]=loglikelihood(type,data,rates,params,sigma,errtype,nondim,divdatapoint,prcalc,prmlocs,fh1lengths)
    % calulcates log likelihood values for input parameters

    % SSE=sum((1.5-params).^2);
    % sigma2=sigma(1)^2;
    % logll=(-SSE/(2*sigma2))-(0.5*length(params)*log(2*pi*sigma2));
    % logll=-SSE;
    % return

    inputparams=10.^(params);

    % have all but rcap_exp be in log space
    if nondim
        inputparams(3)=params(3);
    else
        inputparams(4)=params(4);
    end

    % have all but delivery location be in log space
    if prcalc
        inputparams(end-1:end)=params(end-1:end);
    end
    
    kpolys=calckpolys(type,rates,inputparams,nondim,prcalc,prmlocs,fh1lengths);

    sigma=10.^sigma;
    
    if nondim
        divvalue=kpolys{divdatapoint}(1,2);
        data(divdatapoint,:)=[];
        kpolys(divdatapoint)=[];
        expdata=data.ratiovalues;
    else
        expdata=data.value;
        divvalue=1;
    end

    simdata=zeros(size(expdata));
    for i=1:length(expdata)
        if data.type(i)=="ratio"
            simdata(i)=kpolys{i}(1,3)/kpolys{i}(1,2);
        elseif nondim
            if data.type(i)=="double"
                simdata(i)=kpolys{i}(1,2);
            elseif data.type(i)=="dimer"
                simdata(i)=kpolys{i}(1,3);
            elseif data.type(i)=="single"
                simdata(i)=kpolys{i}(1,1);
            else
                error('Error. \nNo valid experimental data type.')
            end
            simdata(i)=simdata(i)/divvalue;
        else
            if data.type(i)=="double"
                simdata(i)=kpolys{i}(1,2);
            elseif data.type(i)=="dimer"
                simdata(i)=kpolys{i}(1,3);
            elseif data.type(i)=="single"
                simdata(i)=kpolys{i}(1,1);
            else
                error('Error. \nNo valid experimental data type.')
            end
        end
    end
    logll=0;
    if errtype==1
        types=unique(data.type);
        if ~length(types)==length(sigma)
            error("Number of unique data types not the same as number of sigmas")
        end
        for i=1:length(sigma)
            rows = data.type == types(i);
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

function kpolys=calckpolys(type,rates,params,nondim,prcalc,prmlocs,fh1lengths)
    % calulcates kpolys for input parameters
    if prcalc
        x=params(end-1);
        y=params(end);
        prdobs=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"double"),prmlocs,fh1lengths,'UniformOutput',false);
        prdims=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"dimer"),prmlocs,fh1lengths,'UniformOutput',false);
        for i=1:length(rates.k_delbase)
            prdob=prdobs{i};
            prdim=prdims{i};
            vals=rates.k_delbase{i};
            for j=1:size(vals,1)
                vals(j,1)=vals(j,1)*prdob(j);
                vals(j,2)=vals(j,2)*prdob(j);
                vals(j,3)=vals(j,3)*prdob(j);
                vals(j,4)=vals(j,4)*prdim(j);
                vals(j,5)=vals(j,5)*prdim(j);
            end
            rates.k_delbase{i}=vals;
        end
    end

    % calculate per PRM rates
    if nondim
        % params = alpha_del, deta_cap, rcapp_exp, (gamma_del, tau_rel)
        kcaps=cellfun(@(x) x, rates.k_capbase,'UniformOutput',false); 
        kdels=cellfun(@(x) x.*params(1), rates.k_delbase,'UniformOutput',false); 
        rcaps=cellfun(@(x) ((x).^params(3)).*params(2), rates.r_capbase,'UniformOutput',false); 
    else
        kcaps=cellfun(@(x) x.*params(1), rates.k_capbase,'UniformOutput',false);
        kdels=cellfun(@(x) x.*params(2), rates.k_delbase,'UniformOutput',false);
        rcaps=cellfun(@(x) ((x).^params(4)).*params(3), rates.r_capbase,'UniformOutput',false);
    end
    if type=="4st"
        if nondim
            rdels=cellfun(@(x) x.*params(4), rates.r_delbase,'UniformOutput',false);
            krels=cellfun(@(x) x.*params(5), rates.k_relbase,'UniformOutput',false);
        else
            rdels=cellfun(@(x) x.*params(5), rates.r_delbase,'UniformOutput',false);
            krels=cellfun(@(x) x.*params(6), rates.k_relbase,'UniformOutput',false);
        end
        kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
    elseif type=="3st"
        rdels=kcaps;
        krels=kcaps;
        kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
    end
    

    for i=1:length(kpolys)
        PRMsum=sum(kpolys{i},1); % sum up PRMs
        kpolys{i}=[PRMsum(1),sum(PRMsum(2:3)),sum(PRMsum(4:5))]; % Sum up filaments
    end
end

function trueparams=gettrueparams(params,alphakp,kp,nparams,prcalc)
    %must be 3 state method
    kcap=kp/alphakp;
    trueparams=log10((10.^params)*kcap);
    trueparams(3)=params(3); %rcap_exp
    if prcalc
        trueparams(nparams-1:nparams)=params(nparams-1:nparams); %delivery locations
    end

    i=nparams;
    while i<length(trueparams)
        % dont change sigma values
        i=i+1;
        trueparams(i)=params(i);
    end

    trueparams=[log10(kcap), trueparams];

end