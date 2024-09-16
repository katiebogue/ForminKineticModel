function fig=PRMplot(forminList,parameter,xlab,scale,settings,save)
% PRMplot  Creates 4 scatterplots of polymerization vs. a parameter for
% each PRM.
    %
    %   fig = PRMPLOT(formins,x,'xlab','scale',settings) creates 4
    %   scatterplots of polymerization (scaled as specified) vs. x
    %
    %   fig = PRMPLOT(formins,x,'xlab','scale',settings,true) creates and 
    %   saves 4 scatterplots of polymerization (scaled as specified) vs. x
    %
    %   if saving, it is recomended to set groot 'defaultfigureposition' to 
    %   [400 250 900 750] in order to avoid the figure being cut off when 
    %   saving as a pdf.
    %   
    %   Creates 4 scatterplots:
    %       1. All PRM data, labeled by submodel
    %       2. Single submodel data, labeled by formin
    %       3. Double submodel data, labeled by formin
    %       4. Dimer submodel data, labeled by formin
    %
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization (must be a property in the PRM
    %                     class)
    %       xlab        : x-axis (parameter) label (string)
    %       scale       : kpoly axis scale (can be none, log2, log10, ln)
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %   
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, PRM, OPTIONS, EXPERIMENT, KPOLYPLOT, KPOLYMERIZATION.
    arguments
        forminList Formin
        parameter string
        xlab string
        scale string {mustBeMember(scale,{'none','log2','log10','ln'})} % scale for kpoly
        settings Options
        save logical=false
    end

    if scale=="none"
        yscale=@(x) x;
        ylab="k_{poly}";
    elseif scale=="log2"
        ylab="log_{2}(k_{poly})";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}(k_{poly})";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln(k_{poly})";
        yscale=@(x) log(x);
    end

    title_=append('Polymerization Rates vs. ',xlab,' per individual PRM');
    %% Plot all together
    fig(1)=figure;
    PRMList=[forminList.PRMList];
    kpoly=[PRMList.kpoly];
    kpolydouble=[kpoly.double];
    kpolydimer=[kpoly.dimer];
    xdata=[PRMList.(parameter)];
    kpoly1_individual_scatter = scatter(xdata,yscale([kpoly.single]), 'filled','o');
    hold on
    kpoly2a_individual_scatter = scatter(xdata,yscale([kpolydouble.a]), 'filled','s');
    hold on
    kpoly2b_individual_scatter = scatter(xdata,yscale([kpolydouble.b]), 'filled','s');
    hold on
    kpoly3a_individual_scatter = scatter(xdata,yscale([kpolydimer.a]), 'filled','p');
    hold on
    kpoly3b_individual_scatter = scatter(xdata,yscale([kpolydimer.b]), 'filled','p');

    xlabel(xlab)
    ylabel(ylab)
    legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

    title(title_)

    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end
    
    
    %% Plot individual submodels
    fig(2)=figure;
    typeplot("single")
    fig(3)=figure;
    typeplot("double")
    fig(4)=figure;
    typeplot("dimer")

    function typeplot(type)

        for i=1:length(forminList)
            x=[forminList(i).PRMList.(parameter)];
            y=[forminList(i).PRMList.kpoly];
            y=[y.(type)];
            if type~="single"
                y=[y.a,y.b];
                x=[x,x];
            end
            kpoly_individual_scatter = scatter(x,yscale(y),'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
            hold on
        end
    
        xlabel(xlab)
        ylabel(ylab)
        
        legend([forminList.name], 'Location','northeastoutside');
        
        title(append(title_,' (',type,')'))
         
        if save
            figuresave(gcf,settings,append(gca().Title.String,'.fig'));
        end
        
    end

end