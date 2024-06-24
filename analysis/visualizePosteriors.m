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

    %m = matfile(fullfile(foldername,filename),'Writable',true);
    m = matfile(fullfile(foldername,filename),'Writable',false);

    initialparams=m.params_all_trun(1,:);

    if m.type=="3st" && length(initialparams)<7
        parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp"];
    else
        parameter_names=["k_{cap}","k_{del}","r_{cap}","r_{cap} exp","r_{del}","k_{rel}"];
    end

    n_params=m.nparams;
    for i=1:m.nsigma
        parameter_names=[parameter_names,strcat("sigma",int2str(i))];
    end

    
    fname='Dotum';	fsize = 8;	lw = 3;
    
    % Parameter histograms
    p = 1; NBINS=30;
    figure('units','centimeters','position',[5,5,25,20],'Name','Parameter histograms');hold on;
    BinEdges=zeros(NBINS+1,n_params);
    for i=1:n_params
        subplot(ceil(sqrt(n_params)),ceil(sqrt(n_params)),p); hold on;
        histplot = histogram(m.parameters_all(:,i),NBINS);
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
    containsInf=m.containsInf;
	    figure('units','centimeters','position',[5,5,35,30],'Name','Particle clouds');hold on;
        tiledlayout(n_params,n_params,'TileSpacing','tight','Padding','none')
	    for i= 1:n_params
		    for j = 1:n_params
			    if i > j
				    nexttile(p); hold on;
                    if containsInf(i,1) || containsInf(j,1) 
                        [Ncount,BinCentre] = hist3([m.parameters_all(:,i) m.parameters_all(:,j)],"Edges",{BinEdges(:,i),BinEdges(:,j)});
                    else
                        [Ncount,BinCentre] = hist3([m.parameters_all(:,i) m.parameters_all(:,j)],[NBINS NBINS]);
                    end
                    %BinCentre;
                    if length(unique(BinCentre{1,1}))==length(Ncount) && length(unique(BinCentre{1,2}))==length(Ncount)
			            [c,h]=contourf(BinCentre{1,2},BinCentre{1,1},Ncount,NCON,'linestyle','none');
                    end
                    plot(modeHist(j),modeHist(i),'*','Color',[140, 0, 186]./255,'MarkerSize',8);
			        colorbar('off');%colorbar;grid off;
			        h1=gca;set(h1,'FontName',fname,'FontSize',fsize);
			        if i ~= n_params
				        set(gca,'XTickLabel',[])
			        end
			        if j ~= 1
				        set(gca,'YTickLabel',[])
			        end
			        if i == n_params
				        xlabel(strcat("log_{10} ", parameter_names(j)),'FontName',fname,'FontSize',fsize);
			        end
			        if j == 1
				        ylabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize)
                    end
			    elseif i == j
				    nexttile(p); hold on;
				    title(parameter_names(i))
                    histogram(m.parameters_all(:,i),NBINS);
                    plot(modeHist(i),0,'*','Color',[140, 0, 186]./255,'MarkerSize',8);
                    
                    h1=gca;set(h1,'FontName',fname,'FontSize',fsize);grid on;
				    if i ==n_params
					    xlabel(strcat("log_{10} ",parameter_names(i)),'FontName',fname,'FontSize',fsize);
				    end
				    if j ==1
					    ylabel('Prob. Density','FontName',fname,'FontSize',fsize);
				    end
				    if i ~= n_params
					    %set(gca,'XTickLabel',[])
				    end
				    if j ~= 1
					    set(gca,'YTickLabel',[])
				    end
				    %xlabel(parameter_names(i),'FontName',fname,'FontSize',fsize);%ylabel('Prob. Density','FontName',fname,'FontSize',fsize);
				    axis tight;
				    %xlim([lb(j) ub(j)])
			    end
			    p = p + 1;
		    end
	    end
	    %legend('Param. distr.','Param. best part.')
        set(gcf,'Color','w');
        if(saveTF)
            %savefig(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled_compact.fig'),'compact');
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.png'),'png');
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterCorrelations_Labeled.eps'),'epsc');
        end

        % Parameter VS first 3*NTCHECK iteration
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs first 3*NTCHECK iterations');hold on;
        for i=1:n_params
            subplot(ceil(sqrt(n_params)),ceil(sqrt(n_params)),p);
            plot(m.logparams_all_trun(1:3*m.NTCHECK,i),'LineWidth',0.05);
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
        for i=1:n_params
            subplot(ceil(sqrt(n_params)),ceil(sqrt(n_params)),p);
            plot(m.parameters_all(:,i),'LineWidth',0.05);
            xlabel('Iteration','FontName',fname,'FontSize',fsize);
            ylabel(strcat("log_{10} ", parameter_names(i)),'FontName',fname,'FontSize',fsize);
            p=p+1;
        end
        if(saveTF)
            %saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.eps'),'epsc');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations_Last.png'),'png');
        end


        % Parameter VS iteration
        p = 1;
        figure('units','centimeters','position',[5,5,25,20],'Name','Parameter vs iteration');hold on;
        
        for i=1:n_params
            subplot(ceil(sqrt(n_params)),ceil(sqrt(n_params)),p);
            %scatter(1:maxrow,m.params_all(1:maxrow,i));
            plot(m.logparams_all_trun(:,i),'LineWidth',0.05);
            xlabel('Iteration','FontName',fname,'FontSize',fsize);
            ylabel(strcat("log_{10}", parameter_names(i)),'FontName',fname,'FontSize',fsize);
            p=p+1;
        end
        if(saveTF)
            %saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.fig'),'fig');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.eps'),'epsc');
            saveas(gcf,fullfile(savefigfolder,'ParameterVSIterations.png'),'png');
        end


        if (plotecdfTF)
            % Parameter ecdf
            p = 1; NBINS=30;
            figure('units','centimeters','position',[5,5,25,20],'Name','Parameter ecdfs');hold on;
            for i=1:n_params
                % calculate ecdf
                [f_ecdf,x_ecdf] = ecdf(m.parameters_all(1:sizeparameters_all,i));
    
                % calculate ci
                ci_low_ind = find(f_ecdf<0.025,1,'last');
                p_ci(i,1) = x_ecdf(ci_low_ind);
    
                ci_up_ind = find(f_ecdf>0.975,1,'first');
                p_ci(i,2) = x_ecdf(ci_up_ind);
    
                % plot, with ci
                subplot(ceil(sqrt(n_params)),ceil(sqrt(n_params)),p); hold on;
    
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
