function fig = forminkpolybar(formin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    formin Formin
end
% creates polymerization bar charts for all three states

    colors=formin.opts.colors;

    fig=figure; 

    rateplot("kcap",1,"k_{cap}");
    rateplot("kdel",2,"k_{del}");
    rateplot("kpoly",3,"k_{poly}");
    
    function rateplot(rate,plotn,title_)
        rates=[formin.PRMList.(rate)];
        ratedouble=[rates.double];
        ratedimer=[rates.dimer];
        rateplot=[[rates.single];[ratedouble.a]+[ratedouble.b];[ratedimer.a]+[ratedimer.b]];
        fil=[1 2 3];
        
        subplot(1,3,plotn) %1x3 grid w/ axis on 1st cell
        
        rate_bar = bar(fil,rateplot,0.5, 'stacked');
        title(title_);
        ylabel('rate (s^{-1})');
        xticklabels({'single', 'double' , 'dimer'});
        xtickangle(45);
        
        for i = 1:length(rates)
            rate_bar(i).FaceColor = colors(i);
        end
    end
end