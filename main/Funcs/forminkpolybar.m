function fig = forminkpolybar(formin,scale)
% FORMINKPOLYBAR  creates bargraphs of the rates of each step in and of 
% polymerization for the formin.
    %   fig = FORMINKPOLYBAR(formin) creates figure with subplot for the 
    %   rate of each step calculated for formin. 
    %
    %   fig = FORMINBAR(formin,'scale') creates figure with subplot for the
    %   rate of each step (scaled as specified) calculated for formin. 
    %
    %   Determines the steps from formin.opts.equations
    %   Bargraph is stacked, colored by PRM
    %
    %   Inputs:
    %       formin : Formin object
    %       scale  : kpoly axis scale (can be none, log2, log10, ln)
    %   
    %   See also FORMIN, PRM, EXPERIMENT, OPTIONS.
arguments
    formin Formin
    scale {mustBeMember(scale,{'none','log2','log10','ln'})}="none"
end
% creates polymerization bar charts for all three states

    colors=formin.opts.colors;

    fig=figure; 

    if scale=="none"
        yscale=@(x) x;
        ylab="rate (s^{-1})";
    elseif scale=="log2"
        ylab="log_{2}(rate (s^{-1}))";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}(rate (s^{-1}))";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln(rate (s^{-1}))";
        yscale=@(x) log(x);
    end

    ratenames=fieldnames(formin.opts.equations);
    numrates=length(ratenames);
    numsubplots=numrates+1;
    for j=1:numrates
        rateplot(ratenames{j},j,strcat(insertAfter(ratenames{j},1,"_{"),"}"))
    end
    rateplot("kpoly",numsubplots,"k_{poly}");
    
    function rateplot(rate,plotn,title_)
        rates=[formin.PRMList.(rate)];
        if class(rates)=="FilType"
            ratedouble=[rates.double];
            ratedimer=[rates.dimer];
            rateplot=[yscale([rates.single]);yscale([ratedouble.a]+[ratedouble.b]);yscale([ratedimer.a]+[ratedimer.b])];
        else
            rateplot=[yscale(rates);yscale(rates);yscale(rates)];
        end
        fil=[1 2 3];
        
        subplot(1,numsubplots,plotn) %1x3 grid w/ axis on 1st cell
        
        rate_bar = bar(fil,rateplot,0.5, 'stacked');
        title(title_);
        if plotn==1
            ylabel(ylab);
        end
        xticklabels({'single', 'double' , 'dimer'});
        xtickangle(45);
        
        for i = 1:length(rates)
            rate_bar(i).FaceColor = colors(i);
        end
    end
end