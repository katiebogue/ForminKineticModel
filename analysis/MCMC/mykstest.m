function [KSstatistic,histplot,cdfplot] = mykstest(x1,x2,x1label,x2label,plotTF,ax)
%MYKSTEST compute Kolmogorov–Smirnov statistic between two distributions 
%
%   [KSstatistic,histplot,cdfplot] = MYKSTEST(x1,x2) compute Kolmogorov–Smirnov statistic between x1
%   and x2 and generate plot with two panels: an overlapping histogram and a cdf plot
%
%   [KSstatistic,histplot,cdfplot] = MYKSTEST(x1,x2,x1label,x2label) compute Kolmogorov–Smirnov statistic between x1
%   and x2, and generate plot with specified x and y labels with two
%   panels: an overlapping histogram and a cdf plot
%
%   [KSstatistic,histplot,cdfplot] = MYKSTEST(x1,x2,x1label,x2label,0) compute Kolmogorov–Smirnov statistic between x1
%   and x2, do not make a plot
%
%   [KSstatistic,histplot,cdfplot] = MYKSTEST(x1,x2,x1label,x2label,1,ax) compute Kolmogorov–Smirnov statistic between x1
%   and x2, and generate plot on axes ax with specified x and y labels with two
%   panels: an overlapping histogram and a cdf plot
%
%   Inputs:
%       x1       : array of values to compare with x2
%       x2       : array of values to compare with x1
%       x1label  : x axes label for x1 values (default is "x1")
%       x2label  : x axes label for x2 values (default is "x2")
%       plotTF   : whether or not to generate plot (default is true)
%       ax       : 2x1 axes object to add plots to (default is to just
%                   generate new figures)
%       
%   Outputs:
%       KSstatistic : KS stat value
%       histplot    : histogram object with generated histograms
%       cdfplot     : cdf plot object
%
% See also MCMCPARAMFIT.
arguments
    x1
    x2
    x1label="x1"
    x2label="x2"
    plotTF=1
    ax=0 %should be a 2x1 axes object
end

    %
    % Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
    %

    [N,binEdges]=histcounts([x1;x2]);

    binCounts1  =  histcounts (x1 , binEdges);
    binCounts2  =  histcounts (x2 , binEdges);
    
    sumCounts1  =  cumsum(binCounts1)./length(x1);
    sumCounts2  =  cumsum(binCounts2)./length(x2);
    
    sampleCDF1  =  sumCounts1(1:end-1);
    sampleCDF2  =  sumCounts2(1:end-1);
    
    %
    % Compute the test statistic of interest.
    %
    
    %  2-sided test: T = max|F1(x) - F2(x)|.
    deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
    
    
    [KSstatistic,index]   =  max(deltaCDF);

    if plotTF
        if isnumeric(ax)
            figure('units','centimeters','position',[5,5,35,20],'Name','KStest');hold on;
            tiles = tiledlayout(1,2,'TileSpacing','tight','Padding','none');
            ax1=nexttile(1);
            ax2=nexttile(2);
        else
            ax1=ax(1);
            ax2=ax(2);
        end

        hist1  =  histogram (ax1,x1 , binEdges);
        hold(ax1,'on')
        hist2  =  histogram (ax1,x2 , binEdges);
        legend(ax1,{x1label,x2label})
        ylabel(ax1,"Counts")
        histplot=ax1;

        plot(ax2,sampleCDF1)
        hold(ax2,'on')
        plot(ax2,sampleCDF2)
        plot(ax2,deltaCDF)
        plot(ax2,index,KSstatistic,'*','MarkerSize',8);

        legend(ax2,{strcat(x1label," eCDF"),strcat(x2label," eCDF"),"deltaCDF",strcat("KS=",num2str(KSstatistic))})
        if KSstatistic<0.01
            xregion(ax2,index-1,index+1,'HandleVisibility','off','FaceColor','green')
        elseif KSstatistic<0.05
            xregion(ax2,index-1,index+1,'HandleVisibility','off','FaceColor','#EDB120')
        else
            xregion(ax2,index-1,index+1,'HandleVisibility','off','FaceColor','red')
        end
        cdfplot=ax2;

        hold(ax1,'off')
        hold(ax2,'off')

    end

end