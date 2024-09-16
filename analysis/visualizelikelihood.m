function visualizelikelihood(foldername,fitparams)
% please make sure fitparams are entered in log form and with sigma
% values and not in nondimensionalized form

    filename='mcmc_results.mat';

    load(fullfile(foldername,filename),'nparams')
    load(fullfile(foldername,filename),'nsigma')
    load(fullfile(foldername,filename),'containsInf')
    load(fullfile(foldername,filename),'NTCHECK')
    load(fullfile(foldername,filename),'type')
    load(fullfile(foldername,filename),'nondim')
    load(fullfile(foldername,filename),'parameters_all')
    load(fullfile(foldername,filename),'minlogll_params')
    load(fullfile(foldername,filename),'nkpolyparams')


    [minlogll_logll,minlogll_SSE]=getlogll(minlogll_params(1:nkpolyparams+nondim),type,10.^minlogll_params(nkpolyparams+1+nondim:end),nondim)
    [fitparams_logll,fitparams_SSE]=getlogll(fitparams(1:nkpolyparams+nondim),type,10.^fitparams(nkpolyparams+1+nondim:end),nondim)

    load('Users/katiebogue/MATLAB/GitHub/kpolyMCMC/Experiments_4c.mat')
    Experiment1.opts.set_equation(1);
    Experiment1.opts.k_cap=10^fitparams(1);
    Experiment1.opts.k_del=10^fitparams(2);
    Experiment1.opts.r_cap=10^fitparams(3);
    Experiment1.opts.r_cap_exp=fitparams(4);
    Experiment1.opts.kpoly_type=type;
    
    if type=="4st"
        Experiment1.opts.r_del=10^fitparams(5);
        Experiment1.opts.k_rel=10^fitparams(6);
    end
    Experiment1.alldataplot;

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
    Experiment1.alldataplot;


    

    if type=="3st" && length(parameters_all(1,:))<7
        if nondim
            parameter_names=["alpha_{del}","beta_{cap}","r_{cap} exp"];
            fitparamstemp=log10((10.^(fitparams))./10^fitparams(1));
            fitparams=[fitparamstemp(2),fitparamstemp(3),fitparams(4),fitparams(5:end)];
        else
            parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp"];
        end
    else
        if nondim
            parameter_names=["alpha_{del}","beta_{cap}","r_{cap} exp","gamma_{del}","tau_{rel}"];
            fitparamstemp=log10((10.^(fitparams))./10^fitparams(1));
            fitparams=[fitparamstemp(2),fitparamstemp(3),fitparams(4),fitparamstemp(5),fitparamstemp(6),fitparams(7:end)];
        else
            parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp","r_{del}","k_{rel}"];
        end
        
    end

    for i=1:nsigma
        parameter_names=[parameter_names,strcat("sigma",int2str(i))];
    end

    parameter_labels=parameter_names;
    for i=1:length(parameter_names)
        if parameter_names(i)~="r_{cap} exp"
            parameter_labels(i)=strcat("log_{10} ", parameter_names(i));
        end
    end
    fname='Dotum';	fsize = 8;	lw = 3;

    p=1; NCON=300; NBINS=30;
    BinEdges=zeros(NBINS+1,nparams);
    for i=1:nparams
        [N,edges] = histcounts(parameters_all(:,i),NBINS);

        % Find modes of histograms
        BinEdges(:,i)=edges;
        [maxVal, maxInd] = max(N);
        modeHist(i) = edges(maxInd)+0.5*(edges(2)-edges(1));
    end

    % Parameter correlations
    % 2D plots

    if nondim
        load(fullfile(foldername,filename),'minlogll_params_raw')
        minlogll_params=minlogll_params_raw;
    end

    figure('units','centimeters','position',[5,5,35,30],'Name','Particle clouds');hold on;
    tiledlayout(nparams,nparams,'TileSpacing','tight','Padding','none')
    for i= 1:nparams
	    for j = 1:nparams
		    if i > j
			    nexttile(p); hold on;
                if containsInf(i,1) || containsInf(j,1) 
                    [Ncount,BinCentre] = hist3([parameters_all(:,i) parameters_all(:,j)],"Edges",{BinEdges(:,i),BinEdges(:,j)});
                else
                    [Ncount,BinCentre] = hist3([parameters_all(:,i) parameters_all(:,j)],[NBINS NBINS]);
                end
                %BinCentre;
                if length(unique(BinCentre{1,1}))==length(Ncount) && length(unique(BinCentre{1,2}))==length(Ncount)
		            [c,h]=contourf(BinCentre{1,2},BinCentre{1,1},Ncount,NCON,'linestyle','none');
                end
                plot(modeHist(j),modeHist(i),'*','Color',[140, 0, 186]./255,'MarkerSize',8);
                plot(fitparams(j),fitparams(i),'o','MarkerFaceColor','red','MarkerSize',9);
                plot(minlogll_params(j),minlogll_params(i),'o','MarkerFaceColor','blue','MarkerSize',9);
		        colorbar('off');%colorbar;grid off;
		        h1=gca;set(h1,'FontName',fname,'FontSize',fsize);
		        if i ~= nparams
			        set(gca,'XTickLabel',[])
		        end
		        if j ~= 1
			        set(gca,'YTickLabel',[])
		        end
		        if i == nparams
			        xlabel(parameter_labels(j),'FontName',fname,'FontSize',fsize);
		        end
		        if j == 1
			        ylabel(parameter_labels(i),'FontName',fname,'FontSize',fsize)
                end
		    elseif i == j
			    nexttile(p); hold on;
			    title(parameter_names(i))
                
                histogram(parameters_all(:,i),NBINS,'HandleVisibility','off');

                plot(modeHist(i),0,'*','Color',[140, 0, 186]./255,'MarkerSize',8);
                plot(fitparams(i),0,'o','MarkerFaceColor','red','MarkerSize',9);
                plot(minlogll_params(i),0,'o','MarkerFaceColor','blue','MarkerSize',9);
                
                h1=gca;set(h1,'FontName',fname,'FontSize',fsize);grid on;
			    if i ==nparams
				    xlabel(parameter_labels(i),'FontName',fname,'FontSize',fsize);
			    end
			    if j ==1
				    ylabel('Prob. Density','FontName',fname,'FontSize',fsize);
			    end
			    if i ~= nparams
				    %set(gca,'XTickLabel',[])
			    end
			    if j ~= 1
				    set(gca,'YTickLabel',[])
			    end
			    axis tight;
		    end
		    p = p + 1;
	    end
    end
    %legend('Param. distr.','Param. best part.')
    lgd=legend("Histogram Mode",strcat("Test Params, logll=",num2str(fitparams_logll)),strcat("MCMC maxll Params, logll=", num2str(minlogll_logll)));
    lgd.Layout.Tile = 2;
    set(gcf,'Color','w');