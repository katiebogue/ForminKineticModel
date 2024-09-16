function MCMC_pred(Exp,mcmcfile,saveTF)
% 
    %   out = 
    %   
    %   Inputs:
    %         : 
    %
    %   
    %   
    %   See 


    NBINS=30;

    opts=Exp.opts;
    out_struct=readinExp(Exp);
    rates=out_struct.rates;
    clear Exp

    params=10.^(trueparamconvert(mcmcfile));

    load(mcmcfile,'type')

    kpoly_single=zeros(length(params),1);
    kpoly_double=zeros(length(params),1);
    kpoly_dimer=zeros(length(params),1);

    for i=1:length(params)
        kpolys=calckpolys(type,rates,params(i,:));
        kpoly_single(i,:)=kpolys(1);
        kpoly_double(i,:)=kpolys(2);
        kpoly_dimer(i,:)=kpolys(3);
    end

    kpoly_ratios=kpoly_dimer./kpoly_double;

    figure('units','centimeters','position',[5,5,30,25],'Name','Kpoly histograms');hold on;
    %title("BNI1 K_{poly} predictions from MCMC, 3 state")
    histsingle= makehistplot(kpoly_single,NBINS,"K_{poly} Single",1);
    histdouble= makehistplot(kpoly_double,NBINS,"K_{poly} Double",2);
    histdimer= makehistplot(kpoly_dimer,NBINS,"K_{poly} Dimer",3);
    histratio= makehistplot(kpoly_ratios,NBINS,"K_{poly} Ratio Dimer/Double",4);

    if(saveTF)
        saveas(gcf,'KpolyHistograms.eps','epsc');
        saveas(gcf,'KpolyHistograms.png','png');
    end

end

function histplot=makehistplot(x,nbins,xlab,loc)
subplot(2,2,loc); hold on;
fname='Dotum';	fsize = 8;	lw = 3;
histplot=histogram(x,nbins,'HandleVisibility','off');
stdev=std(x);
meanx=mean(x);

xlabel(xlab,'FontName',fname,'FontSize',fsize);
ylabel('Frequency','FontName',fname,'FontSize',fsize);
    
    
% Find modes of histograms
[maxVal, maxInd] = max(histplot.Values);
modeHist = histplot.BinEdges(maxInd)+0.5*histplot.BinWidth;

xline(modeHist, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--')
xline(meanx, 'Color', 'green', 'LineWidth', 2, 'LineStyle', '--')

xline(modeHist - stdev, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(modeHist + stdev, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--','HandleVisibility','off');

legend(sprintf('mode= %3.2f', modeHist),sprintf('mean = %3.2f', meanx),sprintf('std = %3.2f',stdev),'Location','northwest') 

end

function kpolys=calckpolys(type,rates,params)
    % calulcates kpolys for input parameters

    % calculate per PRM rates
    kcaps=cellfun(@(x) x.*params(1), rates.k_capbase,'UniformOutput',false);
    kdels=cellfun(@(x) x.*params(2), rates.k_delbase,'UniformOutput',false);
    rcaps=cellfun(@(x) ((x).^params(4)).*params(3), rates.r_capbase,'UniformOutput',false);
    if type=="4st"
        rdels=cellfun(@(x) x.*params(5), rates.r_delbase,'UniformOutput',false);
        krels=cellfun(@(x) x.*params(6), rates.k_relbase,'UniformOutput',false);
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

    kpolys=kpolys{1};
end

