function fig=PRMplot(forminList,parameter,xlab,settings,save)
% PRMplot  Creates and saves 4 scatterplots of polymerization vs. a
% parameter for each PRM.
    %
    %   PRMPLOT(x,'xlab',Data,settings) creates and saves 4
    %   scatterplots of polymerization vs. x for each PRM
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
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       d           : structure containing data
    %                     must include: 
    %                     all_log_kp1, all_log_kp2a, all_log_kp2b,
    %                     all_log_kp3a, all_log_kp3b, all_iSite_tot,
    %                     all_fh1_names_nobind 
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %   
    %   Calls figuresave
    %
    %   See also .
    arguments
        forminList Formin
        parameter string
        xlab string
        settings Options
        save logical=false
    end

    if save
        set(groot,'defaultfigureposition',[400 250 900 750])
    end

    title_=append('Polymerization Rates vs. ',xlab,' per individual PRM');
    %% Plot all together
    fig(1)=figure;
    PRMList=[forminList.PRMList];
    kpoly=[PRMList.kpoly];
    kpolydouble=[kpoly.double];
    kpolydimer=[kpoly.dimer];
    xdata=[PRMList.(parameter)];
    kpoly1_individual_scatter = scatter(xdata,[kpoly.single], 'filled','o');
    hold on
    kpoly2a_individual_scatter = scatter(xdata,[kpolydouble.a], 'filled','s');
    hold on
    kpoly2b_individual_scatter = scatter(xdata,[kpolydouble.b], 'filled','s');
    hold on
    kpoly3a_individual_scatter = scatter(xdata,[kpolydimer.a], 'filled','p');
    hold on
    kpoly3b_individual_scatter = scatter(xdata,[kpolydimer.b], 'filled','p');

    xlabel(xlab)
    ylabel('(kpoly)')
    legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

    title(title_)
    
    
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
            kpoly_individual_scatter = scatter(x,y,'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
            hold on
        end
    
        xlabel(xlab)
        ylabel('(kpoly)')
        
        legend([forminList.name], 'Location','northeastoutside');
        
        title(append(title_,' (',type,')'))
         
        if save
            figuresave(gcf,settings,append(gca().Title.String,'.fig'));
        end
        
    end

end