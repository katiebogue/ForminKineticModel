function [minlogll_params, minlogll_params_raw,fitTF] = SimulateAnneal(Exp,initialguess,exptype,type,errtype, inputtype, nondim, NTCHECK,NTADAPT,NTMAX,KSCRITICAL)
% 
    %   out = SimulateAnneal(Exp) 
    %   
    %   Inputs:
    %         : 
    %
    %   
    %   
    %   See 

    arguments
        Exp
        initialguess
        exptype % 1= is an experiment, 2= is a struct with rates, data, and resultsfolder and resultsdir
        type
        errtype % 1= use separate sigma values, 2= divide each by SEM
        inputtype % 1= nondimensionalize params, 2= regular params
        nondim= 1 % whether to use nondimensionality
        NTCHECK = 1000
        NTADAPT =100
        NTMAX =10^6
        KSCRITICAL =0.01
    end

    
    NBINS = 200;
    PARAMMAX = 30; % in log-space
    PARAMMIN = -10; % in log-space
    SIGMAMAX = 2; % in log-space
    SIGMAMIN = -2; % in log-space
    EXPMIN = 0.1; % non log-space, applies to the 4th parameter (or 3rd if nondimensional)
    EXPMAX = 2; % non log-space, applies to the 4th parameter (or 3rd if nondimensional)


    % set up place to store data
    if exptype==1
        opts=Exp.opts;
        out_struct=readinExp(Exp);
        rates=out_struct.rates;
        datatab=out_struct.data;
        opts.update_results_folder
        opts.resultsfolder=strcat("MCMC_",opts.resultsfolder);
        clear Exp
        %disp("Read in information from Experiment object")
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
    nparams=nkpolyparams+nsigma;

    % Pick initial guess
    % note: all in log10 space
    rng("shuffle")
    % params=rand([1 nparams])+1;
    % dx=ones(nparams,1);
    if nondim
        if inputtype==2
            initialguess(2:3)=log10(10.^(initialguess(2:3))./10^(initialguess(1)));
            if type=="4st"
                initialguess(5:6)=log10(10.^(initialguess(5:6))./10^(initialguess(1)));
            end
            initialguess(1)=[];
        end
    else
        if inputtype==1
            error("must supply regular parameters")
        end
    end
    params=initialguess;
    dx=ones(nparams,1);
    dx(1:3)=[0.1,0.1,0.1];
    if nondim
        dx(3)=1;
    else
        dx(4)=1;
    end
    if type=="4st"
        if nondim
            dx(5:6)=[1,5];
        else
            dx(5:6)=[5,5];
        end
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

    paccept_matrix(ceil(NTMAX/NTCHECK)+1,nparams+1)=0;
    paramHistCounts_matrices(nparams,NBINS,ceil(NTMAX/NTCHECK))=0;
    paramHistCounts_matrices_nts(ceil(NTMAX/NTCHECK),1)=0;
    ksvals(ceil(NTMAX/NTCHECK),nparams)=0;

    params_temp=zeros(NTCHECK,nparams);
   
    nt_temp=0;
    last_nt=0;
    ntcheck_count=0;
    ksvalindex=1;
    currentntcheck=NTCHECK;
    paramHistCounts = zeros(nparams,NBINS);
    paramHistCountsPrevious = zeros(nparams,NBINS);
    [logll_nt,divvalue]=loglikelihood(type,data,rates,params(1:nkpolyparams),params(nkpolyparams+1:nparams),errtype,nondim,divdatapoint);
    minlogll=logll_nt;
    if nondim
        minlogll_params=gettrueparams(params,divvalue,divkpoly,nkpolyparams);
    else
        minlogll_params=params;
    end
    minlogll_params_raw=params;
    %fitTF=checkfit(type,data,rates,params(1:nkpolyparams),nondim,divdatapoint);

    
    %disp("starting MCMC loop")
    while(nt<NTMAX)
        nt=nt+1;
        nt_temp=nt_temp+1;
        
        if (nt<NTCHECK-NTADAPT)&&(rem(nt,NTADAPT)==0)

            % disp("adapting step size...")
            % fprintf('nt: %d\n',nt)

            proposals_new=proposals-proposals_adapt;
            accepts_new=accepts-accepts_adapt;
            
            p_accept=accepts./proposals;
            p_accept=accepts_new./proposals_new;
            % disp("new acceptance probability:")
            % disp(p_accept')
            p_accept(p_accept==0)=0.01;
            p_accept(proposals==0)=0.44;
            in=(p_accept>0.6 | p_accept<0.3);
            % disp("previous step sizes: ")
            % disp(dx')
            dx(in)=dx(in).*(p_accept(in)./0.44);

            % disp("new step sizes: ")
            % disp(dx')
            proposals_adapt=proposals;
            accepts_adapt=accepts;
        end

        %randomly perturb one parameter
        [proposal,index]=generateproposal(params,dx);
        proposals(index)=proposals(index)+1;
        proposals_temp(index)=proposals_temp(index)+1;

        if (index~=(4-nondim) && proposal(index)<PARAMMIN) || (index~=(4-nondim) && proposal(index)>PARAMMAX) || (index>nkpolyparams && (proposal(index)>SIGMAMAX || proposal(index)<SIGMAMIN)) || (index==(4-nondim) && (proposal(index)>EXPMAX || proposal(index)<EXPMIN))
            % reject anything beyond the boundaries
        else
            %calculate logll of proposal
            [logll_prop,divvalue]=loglikelihood(type,data,rates,proposal(1:nkpolyparams),proposal(nkpolyparams+1:nparams),errtype,nondim,divdatapoint);
            % trueparams=gettrueparams(proposal,divvalue,divkpoly,nkpolyparams);
            % exptrueparams=10.^(trueparams);
            % exptrueparams(4)=trueparams(4);
            % kps=calckpolys(type,rates,exptrueparams,0);
            % val=kps{divdatapoint}(1,2);
            % if abs(val-divkpoly)>10^-3
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
                        minlogll_params=gettrueparams(params,divvalue,divkpoly,nkpolyparams);
                        minlogll_params_raw=params;
                    else
                        minlogll_params=params;
                    end
                end
            elseif rand < exp((logll_prop-logll_nt)/kbt(nt))
                % Boltzmann test, Accept 
                params=proposal;
                logll_nt=logll_prop;
                %kpolys_nt=kpolys_prop;
                accepts(index)=accepts(index)+1;
                accepts_temp(index)=accepts_temp(index)+1;
            end
        end

        if nt<NTCHECK
            params_temp(nt,:)=params;
        end
         if nt==NTCHECK
            binEdges=zeros(nparams,NBINS+1);
            Y = zeros(nparams,NBINS);
            for i=1:nparams
                [Y(i,:),binEdges(i,:)]=histcounts(params_temp((NTCHECK/2):end,i),NBINS);
            end
            paramHistCounts_matrices(:,:,1)=Y;
            paramHistCounts_matrices_nts(1,1)=nt;
            clear Y
            HistCountIndex=2;
        end

        if nt_temp==3*NTCHECK
            ntcheck_count=ntcheck_count+1;
            paccept_matrix(ntcheck_count,1)=nt;
            paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            %disp([accepts_temp./proposals_temp]')
            accepts_temp(:)=0;
            proposals_temp(:)=0;
            last_nt=nt;
            nt_temp=0;
        end

        if nt==NTMAX
            ntcheck_count=ntcheck_count+1;
            paccept_matrix(ntcheck_count,1)=nt;
            paccept_matrix(ntcheck_count,2:end)=[accepts_temp./proposals_temp]';
            %disp([accepts_temp./proposals_temp]')
            accepts_temp(:)=0;
            proposals_temp(:)=0;
            last_nt=nt;
            nt_temp=0;
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
        end

        if nt==3*currentntcheck
            %fprintf('nt: %d\n',nt)
            %fprintf('currentntcheck: %d\n',currentntcheck)

            paramHistCounts_matrices_nts(HistCountIndex,1)=currentntcheck*2;
            paramHistCounts_matrices_nts(HistCountIndex+1,1)=currentntcheck*3;
            paramHistCounts_matrices(:,:,HistCountIndex)=paramHistCountsPrevious;
            paramHistCounts_matrices(:,:,HistCountIndex+1)=paramHistCounts;
            HistCountIndex=HistCountIndex+2;
            
            %disp('Begining KS test')
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
            ksvals(ksvalindex,:)=ksStatistic';
            ksvalindex=ksvalindex+1;
            if all(ksStatistic(:) < KSCRITICAL)
                %disp('KS test successful')
                maxlikelihoodplots(minlogll_params,type)
                fitTF=checkfit(type,data,rates,params(1:nkpolyparams),nondim,divdatapoint);
                return
            else
                %disp(ksStatistic')
                %disp('KS test unsuccessful')
                currentntcheck=currentntcheck*3;
                %fprintf('new currentntcheck: %d\n',currentntcheck)
                paramHistCounts = zeros(nparams,NBINS);
                paramHistCountsPrevious = zeros(nparams,NBINS);
            end
        end
    end

    maxlikelihoodplots(minlogll_params,type)
    fitTF=checkfit(type,data,rates,params(1:nkpolyparams),nondim,divdatapoint);

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

    function output = kbt(ntvalue)
        %output=1-((ntvalue/NTMAX)^(8));

        output=0.5*(NTMAX^5)*(1+cos(ntvalue*pi/NTMAX));
    end

end

function maxlikelihoodplots(minlogll_params,type)
    load('Users/katiebogue/MATLAB/GitHub/kpolyMCMC/Experiments_4c.mat')

    Experiment1.opts.set_equation(1);
    Experiment1.opts.k_cap=10^minlogll_params(1);
    Experiment1.opts.k_del=10^minlogll_params(2);
    Experiment1.opts.r_cap=10^minlogll_params(3);
    Experiment1.opts.r_cap_exp=minlogll_params(4);
    Experiment1.opts.kpoly_type=type;
    
    if type=="4st"
        Experiment1.opts.r_del=10^minlogll_params(5);
        Experiment1.opts.k_rel=10^minlogll_params(6);
    end
    
    figure('units','centimeters','position',[5,5,45,30],'Name','BestFit');hold on;
    tiles = tiledlayout(2,3,'TileSpacing','tight','Padding','none');
    title(tiles,num2str(minlogll_params))

    fig3=Experiment1.expdatabar(group="Fig 3");
    ax_temp = copyobj(copyobj(fig3(1).Children, get(fig3(1).Children,'Type')), tiles);
    ax_temp(2).Layout.Tile=5;
    close(fig3)

    
    fig35=Experiment1.expdatabar(group="Fig 3 5");
    ax_temp = copyobj(copyobj(fig35(1).Children, get(fig35(1).Children,'Type')), tiles);
    ax_temp(2).Layout.Tile=4;
    close(fig35)

    fig4a=Experiment1.expdatabar(group="Fig 4a");
    ax_temp = copyobj(copyobj(fig4a(1).Children, get(fig4a(1).Children,'Type')), tiles);
    ax_temp(2).Layout.Tile=3;
    close(fig4a)

    fig4c=Experiment1.expdatabar(group="Fig 4c");
    ax_temp = copyobj(copyobj(fig4c(1).Children, get(fig4c(1).Children,'Type')), tiles);
    ax_temp(2).Layout.Tile=2;
    close(fig4c)

    figNTD=Experiment1.expdatabar("log2",false,group="NTD data"); % uses log2 scale
    ax_temp = copyobj(copyobj(figNTD(1).Children, get(figNTD(1).Children,'Type')), tiles);
    ax_temp(2).Layout.Tile=1;
    close(figNTD)
end


function [logll,divvalue]=loglikelihood(type,data,rates,params,sigma,errtype,nondim,divdatapoint)
    % calulcates energy values for input parameters

    inputparams=10.^(params);

    % have all but rcap_exp be in log space
    if nondim
        inputparams(3)=params(3);
    else
        inputparams(4)=params(4);
    end
    
    kpolys=calckpolys(type,rates,inputparams,nondim);

    sigma=10.^sigma;
    
    if nondim
        divvalue=kpolys{divdatapoint}(1,2);
        data(divdatapoint,:)=[];
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
            if i==divdatapoint
            elseif data.type(i)=="double"
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

            %extra penalty for double data
            % if data.type=="double"
            %     logll=logll+2*logll_i;
            % end

            %extra penalty for NTD data
            % if data.type=="ratio"
            %     logll=logll-SSE;
            % end
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

    % extra penalty for more trends
    % rows = data.type == "double";
    % logll=logll-sum(10.*(abs(diff(simdata(rows))-diff(expdata(rows))))); 

    % extra penalty for fig 4a data
    rows=[];
    for i=1:length(expdata)
        if data.groups{i}=="Fig 4a"
            rows=[rows, i];
        end
    end
    logll=logll-sum(abs(diff(simdata(rows))-diff(expdata(rows))));
    
    if sign(diff(simdata(rows))) ~= sign(diff(expdata(rows)))
        logll=logll-sum((abs(diff(simdata(rows))-diff(expdata(rows)))));
    end


end

function kpolys=calckpolys(type,rates,params,nondim)
    % calulcates kpolys for input parameters

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

function trueparams=gettrueparams(params,alphakp,kp,nparams)
    %must be 3 state method
    kcap=kp/alphakp;
    trueparams=log10((10.^params)*kcap);
    trueparams(3)=params(3); %rcap_exp

    i=nparams;
    while i<length(trueparams)
        % dont change sigma values
        i=i+1;
        trueparams(i)=params(i);
    end

    trueparams=[log10(kcap), trueparams];

end

function fitTF=checkfit(type,data,rates,params,nondim,divdatapoint)
    % calulcates energy values for input parameters

    inputparams=10.^(params);

    % have all but rcap_exp be in log space
    if nondim
        inputparams(3)=params(3);
    else
        inputparams(4)=params(4);
    end
    
    kpolys=calckpolys(type,rates,inputparams,nondim);
    
    if nondim
        divvalue=kpolys{divdatapoint}(1,2);
        data(divdatapoint,:)=[];
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
            if i==divdatapoint
            elseif data.type(i)=="double"
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

    fitTF=1;

    % check if trend is present for fig 4a
    rows=[];
    for i=1:length(expdata)
        if data.groups{i}=="Fig 4a"
            rows=[rows, i];
        end
    end
    %logll=logll-sum(abs(diff(simdata(rows))-diff(expdata(rows))));
    if sign(diff(simdata(rows))) ~= sign(diff(expdata(rows)))
        fitTF=0;
    end


    % check if trend is present for fig 3
    rows=[];
    for i=1:length(expdata)
        if data.groups{i}=="Fig 3 5"
            rows=[rows, i];
        end
    end
    %logll=logll-sum(abs(diff(simdata(rows))-diff(expdata(rows))));
    if sign(diff(simdata(rows))) ~= sign(diff(expdata(rows)))
        fitTF=0;
    end

    % check if NTD is within error bars
    rows=[];
    for i=1:length(expdata)
        if data.groups{i}=="NTD data"
            rows=[rows, i];
        end
    end
    topvals=data.value(rows)+data.errtop(rows);
    botvals=data.value(rows)-data.errbot(rows);
    for i=rows
        if simdata(i)>topvals(i) || simdata(i)+.01<botvals(i)
            fitTF=0;
        end
    end

end