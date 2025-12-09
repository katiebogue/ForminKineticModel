function kpoly_ratios= MCMC_pred(Exp,mcmcfile,saveTF,boundmin,boundmax)
%MCMC_PRED generate histograms of calculated kpoly values for formin in
%experiment using the parameter cloud from the mcmc file
%
%   MCMC_pred(Exp,mcmcfile) generate histograms of calculated kpoly
%   values for formin in Exp using the parameter cloud in mcmcfile
%
%   MCMC_pred(Exp,mcmcfile,1) generate and save histograms of
%   calculated kpoly values for formin in Exp using the parameter cloud in
%   mcmcfile 
%
%   Inputs:
%       Exp : Experiment file, should only have one formin
%       mcmcfile: .mat file result from mcmcparamfit, must have the
%           following parameters: nparams, nkpolyparams, nsigma, type, nondim,
%           parameters_all, rates, divdatapoint, divkpoly
%       saveTF : whether or not to save the histograms
%
%   Generates a histogram for Kpoly single, double, dimer, and ratio, and
%   includes lines indicating the mean, mode, and standard deviation.
%
%   Runs trueparamconvert.
%
% See also TRUEPARAMCONVERT.
    NBINS=30;

    load(mcmcfile,'prcalc')
    opts=Exp.opts;
    if prcalc
        opts.set_equation(2);
    end
    out_struct=readinExp(Exp);
    rates=out_struct.rates;
    fh1lengths=out_struct.fh1sizes;
    prmlocs=out_struct.prmlocs;
    clear Exp


    if isempty(who('-file',mcmcfile,'vals_trueparamconvert'))
        vals_trueparamconvert = trueparamconvert(mcmcfile);
        save(mcmcfile, 'vals_trueparamconvert',"-append");
    else
        load(mcmcfile, 'vals_trueparamconvert');
    end

    
    params=10.^(vals_trueparamconvert);
    params(4)=vals_trueparamconvert(4); %rcap_exp

    load(mcmcfile,'type')
    load(mcmcfile,'prfit')
    load(mcmcfile,'xloc')
    load(mcmcfile,'yloc')

    kpoly_single=zeros(length(params),1);
    kpoly_double=zeros(length(params),1);
    kpoly_dimer=zeros(length(params),1);

    if prcalc
        if prfit
        else
            x= xloc;
            y= yloc;
            prdobs=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"double",1),prmlocs,fh1lengths,'UniformOutput',false);
            prdims=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"dimer",1),prmlocs,fh1lengths,'UniformOutput',false);
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
    end

    for i=1:length(params)
        kpolys=calckpolys(type,rates,params(i,:),prfit,prmlocs,fh1lengths, prcalc, xloc, yloc);
        kpoly_single(i,:)=kpolys(1);
        kpoly_double(i,:)=kpolys(2);
        kpoly_dimer(i,:)=kpolys(3);
    end

    kpoly_ratios=kpoly_dimer./kpoly_double;

    figure('units','centimeters','position',[5,5,40,25],'Name','Kpoly histograms');hold on;
    t=tiledlayout(2,3,'TileSpacing','tight','Padding','none');
    %title("BNI1 K_{poly} predictions from MCMC, 3 state")
    %histsingle= makehistplot(kpoly_single,NBINS,"K_{poly} Single",1);
    histdouble_cut= makehistplot(kpoly_double,NBINS,"K_{poly} Double",1,1);
    histdimer_cut= makehistplot(kpoly_dimer,NBINS,"K_{poly} Dimer",2,1);
    histratio= makehistplot(kpoly_ratios,NBINS,"K_{poly} Ratio Dimer/Double",3,0);
    if boundmax>histratio.BinLimits(1) && boundmin<histratio.BinLimits(2)
        if boundmin<histratio.BinLimits(1)
            boundmin=histratio.BinLimits(1);
        end
        if boundmax>histratio.BinLimits(2)
            boundmax=histratio.BinLimits(2);
        end
        xregion(boundmin,boundmax,'DisplayName',"target")
    end
    histdouble= makehistplot(kpoly_double,NBINS,"K_{poly} Double",4,0);
    histdimer= makehistplot(kpoly_dimer,NBINS,"K_{poly} Dimer",5,0);

    fontsize(16,"points")

    if(saveTF)
        saveas(gcf,'KpolyHistograms.eps','epsc');
        saveas(gcf,'KpolyHistograms.png','png');
    end

end

function histplot=makehistplot(x,nbins,xlab,loc,cutoff)
%subplot(2,2,loc); hold on;
%subplot(2,3,loc); hold on;
nexttile(loc); hold on;
if cutoff
    mask = x(:) >= 30;
    x(mask)=[];
end
fname='Dotum';	fsize = 10;	lw = 3;
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

xline(min(histplot.Data), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');
xline(max(histplot.Data), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-','HandleVisibility','off');

legend(sprintf('mode= %3.2f', modeHist),sprintf('mean = %3.2f', meanx),sprintf('std = %3.2f',stdev),sprintf('limits'),'Location','best') 

end

function kpolys=calckpolys(type,rates,params,prfit,prmlocs,fh1lengths, prcalc, xloc, yloc)
    % calulcates kpolys for input parameters
    if prcalc
        if prfit
            x=params(end-1);
            y=params(end);
            prdobs=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"double",1),prmlocs,fh1lengths,'UniformOutput',false);
            prdims=cellfun(@(n1,fh1length) pr(n1,fh1length,35.5,1,x,y,"dimer",1),prmlocs,fh1lengths,'UniformOutput',false);
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
    end

    % calculate per PRM rates
        %kcaps=cellfun(@(x) x.*params(1), rates.k_capbase,'UniformOutput',false);
        kcaps = cell(size(rates.k_capbase));
        for i = 1:numel(rates.k_capbase)
            kcaps{i} = rates.k_capbase{i} .* params(1);
        end

        %kdels=cellfun(@(x) x.*params(2), rates.k_delbase,'UniformOutput',false);
        kdels = cell(size(rates.k_delbase));
        for i = 1:numel(rates.k_delbase)
            kdels{i} = rates.k_delbase{i} .* params(2);
        end

        %rcaps=cellfun(@(x) ((x).^params(4)).*params(3), rates.r_capbase,'UniformOutput',false);
        rcaps= cell(size(rates.r_capbase));
        for i = 1:numel(rates.r_capbase)
            rcaps{i} = ((rates.r_capbase{i}).^params(4)).*params(3);
        end
    if type=="4st"
            %rdels=cellfun(@(x) x.*params(5), rates.r_delbase,'UniformOutput',false);
            rdels = cell(size(rates.r_delbase));
            for i = 1:numel(rates.r_delbase)
                rdels{i} = rates.r_delbase{i} .* params(5);
            end

            %krels=cellfun(@(x) x.*params(6), rates.k_relbase,'UniformOutput',false);
            krels = cell(size(rates.k_relbase));
            for i = 1:numel(rates.k_relbase)
                kdels{i} = rates.k_relbase{i} .* params(6);
            end

        %kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
        kpolys = cell(size(kcaps));
        for i = 1:numel(kcaps)
            kcap = kcaps{i};
            kdel = kdels{i};
            rcap = rcaps{i};
            rdel = rdels{i};
            krel = krels{i};
            kpolys{i} = 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel)));
        end
    elseif type=="3st"
        % rdels=kcaps;
        % krels=kcaps;
        %kpolys=cellfun(@(kcap,kdel,rcap) 1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap))),kcaps,kdels,rcaps,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
        kpolys = cell(size(kcaps));
        for i = 1:numel(kcaps)
            kcap = kcaps{i};
            kdel = kdels{i};
            rcap = rcaps{i};
            kpolys{i} = 1 ./ ((1 ./ kdel) + ((kdel + rcap) ./ (kdel .* kcap)));
        end
    end
    

    for i=1:length(kpolys)
        PRMsum=sum(kpolys{i},1); % sum up PRMs
        kpolys{i}=[PRMsum(1),sum(PRMsum(2:3)),sum(PRMsum(4:5))]; % Sum up filaments
    end

    kpolys=kpolys{1};
    if kpolys(1)>10000
         return
    end
end

