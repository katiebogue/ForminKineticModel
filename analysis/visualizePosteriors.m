function visualizePosteriors(foldername,saveTF,savefigfolder,filename)
    
    plotecdfTF = 0; % dont make the ecdf plots (causes MATLAB to crash)

    if ~exist('savefigfolder','var')
        mkdir(foldername,"figures")
        savefigfolder=fullfile(foldername,"figures");
    end
    if ~exist('filename','var')
        filename='mcmc_results.mat';
    end

    updateMCMCoutput(fullfile(foldername,filename));
    reloadmatfile(fullfile(foldername,filename))

    load(fullfile(foldername,filename),'nparams')
    load(fullfile(foldername,filename),'nsigma')
    load(fullfile(foldername,filename),'containsInf')
    load(fullfile(foldername,filename),'NTCHECK')
    load(fullfile(foldername,filename),'type')

    inmemTF=[0,0];
    % logparams_all_trun=0;
    % parameters_all=0;
    try
        load(fullfile(foldername,filename),'logparams_all_trun')
        initialparams=logparams_all_trun(1,:);
        inmemTF(1)=1;
    catch
        m = matfile(fullfile(foldername,filename),'Writable',false);
        initialparams=m.logparams_all_trun(1,:);
    end

    try
        load(fullfile(foldername,filename),'parameters_all')
        inmemTF(2)=1;
    catch
        if inmemTF(1)
            m = matfile(fullfile(foldername,filename),'Writable',false);
        end
    end

    % m = matfile(fullfile(foldername,filename),'Writable',false);
    % initialparams=m.logparams_all_trun(1,:);

    if type=="3st" && length(initialparams)<7
        parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp"];
    else
        parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp","r_{del}","k_{rel}"];
    end

    for i=1:nsigma
        parameter_names=[parameter_names,strcat("sigma",int2str(i))];
    end

    
    fname='Dotum';	fsize = 8;	lw = 3;
    
    % Parameter histograms
    p = 1; NBINS=30;
    figure('units','centimeters','position',[5,5,25,20],'Name','Parameter histograms');hold on;
    BinEdges=zeros(NBINS+1,nparams);
    for i=1:nparams
        subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p); hold on;
        if inmemTF(2)
            histplot = histogram(parameters_all(:,i),NBINS);
        else
            histplot = histogram(m.parameters_all(:,i),NBINS);
        end
        xlabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize);
        ylabel('Frequency','FontName',fname,'FontSize',fsize);
    
    
        % Find modes of histograms
        BinEdges(:,i)=histplot.BinEdges;
        [maxVal, maxInd] = max(histplot.Values);
        modeHist(i) = histplot.BinEdges(maxInd)+0.5*histplot.BinWidth;
    
        p=p+1;
    end
    
    if(saveTF)
        %saveas(gcf,fullfile(savefigfolder,'ParameterHistograms.fig'),'fig');
        saveas(gcf,fullfile(savefigfolder,'ParameterHistograms.eps'),'epsc');
        saveas(gcf,fullfile(savefigfolder,'ParameterHistograms.png'),'png');
    end
    
    % Parameter correlations
    % 2D plots
    p=1;	NCON=300;	nX=100;
	    figure('units','centimeters','position',[5,5,35,30],'Name','Particle clouds');hold on;
        tiledlayout(nparams,nparams,'TileSpacing','tight','Padding','none')
	    for i= 1:nparams
		    for j = 1:nparams
			    if i > j
				    nexttile(p); hold on;
                    if containsInf(i,1) || containsInf(j,1) 
                        if inmemTF(2)
                            [Ncount,BinCentre] = hist3([parameters_all(:,i) parameters_all(:,j)],"Edges",{BinEdges(:,i),BinEdges(:,j)});
                        else
                            [Ncount,BinCentre] = hist3([m.parameters_all(:,i) m.parameters_all(:,j)],"Edges",{BinEdges(:,i),BinEdges(:,j)});
                        end
                    else
                        if inmemTF(2)
                            [Ncount,BinCentre] = hist3([parameters_all(:,i) parameters_all(:,j)],[NBINS NBINS]);
                        else
                            [Ncount,BinCentre] = hist3([m.parameters_all(:,i) m.parameters_all(:,j)],[NBINS NBINS]);
                        end
                    end
                    %BinCentre;
                    if length(unique(BinCentre{1,1}))==length(Ncount) && length(unique(BinCentre{1,2}))==length(Ncount)
			            [c,h]=contourf(BinCentre{1,2},BinCentre{1,1},Ncount,NCON,'linestyle','none');
                    end
                    plot(modeHist(j),modeHist(i),'*','Color',[140, 0, 186]./255,'MarkerSize',8);
			        colorbar('off');%colorbar;grid off;
			        h1=gca;set(h1,'FontName',fname,'FontSize',fsize);
			        if i ~= nparams
				        set(gca,'XTickLabel',[])
			        end
			        if j ~= 1
				        set(gca,'YTickLabel',[])
			        end
			        if i == nparams
				        xlabel(strcat("log_{10} ", parameter_names(j)),'FontName',fname,'FontSize',fsize);
			        end
			        if j == 1
				        ylabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize)
                    end
			    elseif i == j
				    nexttile(p); hold on;
				    title(parameter_names(i))
                    if inmemTF(2)
                        histogram(parameters_all(:,i),NBINS);
                    else
                        histogram(m.parameters_all(:,i),NBINS);
                    end

                    plot(modeHist(i),0,'*','Color',[140, 0, 186]./255,'MarkerSize',8);
                    
                    h1=gca;set(h1,'FontName',fname,'FontSize',fsize);grid on;
				    if i ==nparams
					    xlabel(strcat("log_{10} ",parameter_names(i)),'FontName',fname,'FontSize',fsize);
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
        set(gcf,'Color','w');
        if(saveTF)
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.png'),'png');
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.eps'),'epsc');
        end

        % Parameter VS first 3*NTCHECK iteration
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs first 3*NTCHECK iterations');hold on;
        for i=1:nparams
            subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p);
            if inmemTF(1)
                plot(logparams_all_trun(1:3*NTCHECK,i),'LineWidth',0.05);
            else
                plot(m.logparams_all_trun(1:3*NTCHECK,i),'LineWidth',0.05);
            end
            xlabel('Iteration','FontName',fname,'FontSize',fsize);
            ylabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize);
            p=p+1;
        end
        if(saveTF)
            %saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_First.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_First.eps'),'epsc');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_First.png'),'png');
        end

        % Parameter VS posterior iterations
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs posterior iterations');hold on;
        for i=1:nparams
            subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p);
            if inmemTF(2)
                plot(parameters_all(:,i),'LineWidth',0.05);
            else
                plot(m.parameters_all(:,i),'LineWidth',0.05);
            end
            xlabel('Iteration','FontName',fname,'FontSize',fsize);
            ylabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize);
            p=p+1;
        end
        if(saveTF)
            %saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.eps'),'epsc');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.png'),'png');
        end
        close all


        % Parameter VS iteration
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs iteration');hold on;
        
        for i=1:nparams
            subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p);
            if inmemTF(1)
                plot(logparams_all_trun(:,i),'LineWidth',0.05);
            else
                plot(m.logparams_all_trun(:,i),'LineWidth',0.05);
            end
            xlabel('Iteration','FontName',fname,'FontSize',fsize);
            ylabel(strcat("log_{10}", parameter_names(i)),'FontName',fname,'FontSize',fsize);
            p=p+1;
        end
        if(saveTF)
            %saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.eps'),'epsc');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.png'),'png');
        end
        close all


        if (plotecdfTF)
            % Parameter ecdf
            p = 1; NBINS=30;
            figure('units','centimeters','position',[5,5,25,20],'Name','Parameter ecdfs');hold on;
            for i=1:nparams
                % calculate ecdf
                if inmemTF(2)
                    [f_ecdf,x_ecdf] = ecdf(parameters_all(1:sizeparameters_all,i));
                else
                    [f_ecdf,x_ecdf] = ecdf(m.parameters_all(1:sizeparameters_all,i));
                end
    
                % calculate ci
                ci_low_ind = find(f_ecdf<0.025,1,'last');
                p_ci(i,1) = x_ecdf(ci_low_ind);
    
                ci_up_ind = find(f_ecdf>0.975,1,'first');
                p_ci(i,2) = x_ecdf(ci_up_ind);
    
                % plot, with ci
                subplot(ceil(sqrt(nparams)),ceil(sqrt(nparams)),p); hold on;
    
                plot(x_ecdf,f_ecdf,'-','LineWidth',lw);
    
                plot(x_ecdf,0.025*ones(length(x_ecdf),1),'--k');
                plot(x_ecdf,0.975*ones(length(x_ecdf),1),'--k');
    
                plot(p_ci(i,1)*ones(length(f_ecdf),1),f_ecdf,'--r');
                plot(p_ci(i,2)*ones(length(f_ecdf),1),f_ecdf,'--r');
    
                xlabel(parameter_names(i),'FontName',fname,'FontSize',fsize);
                ylabel('Cumulative Distribution','FontName',fname,'FontSize',fsize);
    
                p=p+1;
    
            end
            
            if(saveTF)
                saveas(gcf,fullfile(savefigfolder,'ParameterECDFs.fig'),'fig');
                saveas(gcf,fullfile(savefigfolder,'ParameterECDFs.eps'),'epsc');
            end
        end

        nt=NTCHECK;

        if inmemTF(1)
            len=length(logparams_all_trun);
        else
            len=length(m.logparams_all_trun);
        end

        while nt*3<len
            ksplot(nt)
            nt=nt*3;
        end
        ksplot(floor(len/3))

    function ksplot(nt)
        figure('units','centimeters','position',[5,5,45,30],'Name','KStest');hold on;
        tiles = tiledlayout(nparams,3,'TileSpacing','tight','Padding','none');
    
        q=1;
        ksval=zeros(nparams,1);
        for k=1:nparams
            nexttile(q)
            if inmemTF(1)
                plot(logparams_all_trun(nt:nt*3,k),'LineWidth',0.05);
            else
                plot(m.logparams_all_trun(nt:nt*3,k),'LineWidth',0.05);
            end
            %xlabel('Iteration');
            ylabel(parameter_names(k));
            q=q+1;
    
            ax1=nexttile(q);
            ax2=nexttile(q+1);
            q=q+2;
    
            label1=strcat(num2str(nt),"-",num2str(nt*2));
            label2=strcat(num2str(nt*2+1),"-",num2str(nt*3));
            if inmemTF(1)
                [ksval(k),histploti,cdfplot]=mykstest(logparams_all_trun(nt:nt*2,k),logparams_all_trun(nt*2+1:nt*3,k),label1,label2,1,[ax1,ax2]);
            else
                [ksval(k),histploti,cdfplot]=mykstest(m.logparams_all_trun(nt:nt*2,k),m.logparams_all_trun(nt*2+1:nt*3,k),label1,label2,1,[ax1,ax2]);
            end
        end
    
        tiles.Title.String = strcat("KStest for nt= ",num2str(nt*3)," (ntcheck= ",num2str(nt),")");
    
        if saveTF
            saveas(gcf,fullfile(savefigfolder,strcat("KStestfig_",num2str(nt),".png")'),'png');
            saveas(gcf,fullfile(savefigfolder,strcat("KStestfig_",num2str(nt),".eps")),'epsc');
        end
        close all
    
    end

    load('Users/katiebogue/MATLAB/GitHub/kpolyMCMC/Experiments_4c.mat')
    maxlikelihoodplot(Experiment1,foldername)
        
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions below are called by above function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,xo]=density(x,xout,ss,gaus)
%DENSITY  Density estimator using Gaussian kernel
% Y = DENSITY(X,XOUT,S)
% X is the vector of data values.
% The density estimator is evaluated at XOUT points.
% S is a scale factor for the default kernel bandwidth,
% default S = 1.
% Without output arguments the density is plotted.

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.9 $  $Date: 2012/09/27 11:47:35 $

	if nargin<3
	  ss=1;
	end
	if nargin<4
	  gaus=1;
	end

	if nargin<2 | isempty(xout)
	  xmin=min(x); xmax=max(x); xrange=xmax-xmin;
	  if length(x) > 200
		xout=linspace(xmin-0.08*xrange,xmax+0.08*xrange);
	  else
		xout=linspace(mean(x)-4*std(x),mean(x)+4*std(x));
	  end
	end
	y  = zeros(size(xout));
	n  = length(xout);
	nx = length(x);

	%%% see MASS 2nd ed page 181.
	if iqrange(x)<=0
	  s=1.06*std(x)*nx^(-1/5);
	else
	  s=1.06*min(std(x),iqrange(x)/1.34)*nx^(-1/5);
	end
	%  s=1.144*std(x)*nx^(-1/5);
	if ss>0
	  s=ss*s;
	elseif ss<0
	  s = abs(ss);
	end
	if gaus==1
	  % Gaussian kernel
	  for i=1:n
		y(i) = 1/nx*sum(norpf((xout(i)-x)/s))./s;
%         norpf_in = (xout(i)-x)/s;
%         mu = 0;
%         sigma2 = 1;
%         norpf_out=1./sqrt(2*pi*sigma2).*exp(-0.5*(norpf_in-mu).^2 ./sigma2);
%         y(i) = 1/nx*sum(norpf_out)./s;
      end
	elseif gaus==-1
	  % Gamma kernel (still testing)
	  s = s*0.5;
	  for i=1:n
		y(i) = 1/nx*sum(gammapf(xout(i),x./s+1,s));
	  end
	else
	  % triangular kernel
	  s=s*1.2113;
	  for i=1:n
		y(i) = 1/nx*sum(max(0,1-abs(xout(i)-x)/s))./s;
	  end
	end

	if nargout>1
	  xo=xout;
	end

	if nargout==0
	  plot(xout,y)
	  clear y % no output
	end
end

function y=iqrange(x)
	% Interquantile range of each column of x

	% ML 2000

	% Marko Laine <marko.laine@fmi.fi>
	% $Revision: 1.3 $  $Date: 2012/09/27 11:47:37 $

	[n,m]=size(x);
	if n==1
	  x=x';
	  n = m;
	  m = 1;
	end

	x  = sort(x);
	% n  = length(x);
	i1 = floor((n+1)/4);
	i3 = floor(3/4*(n+1));
	f1 = (n+1)/4-i1;
	f3 = 3/4*(n+1)-i3;
	q1 = (1-f1).*x(i1,:)+f1.*x(i1+1,:);
	q3 = (1-f3).*x(i3,:)+f3.*x(i3+1,:);
	y  = q3-q1;
end

function y=norpf(x,mu,sigma2)
	% NORPF(x,mu,sigma2)  Normal (Gaussian) density function

	% Marko Laine <marko.laine@fmi.fi>
	% $Revision: 1.4 $  $Date: 2012/09/27 11:47:38 $

	if nargin < 2, mu=0; end
	if nargin < 3, sigma2=1; end
	y=1./sqrt(2*pi*sigma2).*exp(-0.5*(x-mu).^2 ./sigma2);
end
