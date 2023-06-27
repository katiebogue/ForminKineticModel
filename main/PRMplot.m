function PRMplot(parameter,xlab,d,settings)
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
    %       parameter   : array containing values to plot against 
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       d           : structure containing data
    %                     must include: 
    %                     all_log_kp1, all_log_kp2a, all_log_kp2b,
    %                     all_log_kp3a, all_log_kp3b, all_iSite_tot,
    %                     all_fh1_names_nobind 
    %       settings    : structure containing settings variables
    %                     must include: 
    %                     settings_variable, pdf_name, shapes, colors
    %   
    %   Calls figuresave
    %
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS, WHOLEPACKAGE, MAKEPOINTS.

title_=append('Polymerization Rates vs. ',xlab,' per individual PRM');
%% Plot all together
kpoly1_individual_scatter = scatter(parameter,d.all_log_kp1, 'filled','o');
hold on
kpoly2a_individual_scatter = scatter(parameter,d.all_log_kp2a, 'filled','s');
hold on
kpoly2b_individual_scatter = scatter(parameter,d.all_log_kp2b, 'filled','s');
hold on
kpoly3a_individual_scatter = scatter(parameter,d.all_log_kp3a, 'filled','p');
hold on
kpoly3b_individual_scatter = scatter(parameter,d.all_log_kp3b, 'filled','p');

xlabel(xlab)
ylabel('log_2(kpoly)')
legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

title(title_)
subtitle(settings.settings_variable)
 
figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));

close all

%% Plot individual submodels
typeplot("single")
typeplot("double")
typeplot("dimer")

    function typeplot(type)
    n=1;
    while n<= length(d.all_log_kp1)
        for LOOPY = 1:length(d.all_iSite_tot)
            PRMnum = d.all_iSite_tot(LOOPY);
            endpoint = n + PRMnum-1;
            if type=="single"
                x=parameter(n:endpoint);
                y=d.all_log_kp1(n:endpoint);
            elseif type=="double"
                x=[parameter(n:endpoint);parameter(n:endpoint)];
                y=[d.all_log_kp2a(n:endpoint);d.all_log_kp2b(n:endpoint)];
            elseif type=="dimer"
                x=[parameter(n:endpoint);parameter(n:endpoint)];
                y=[d.all_log_kp3a(n:endpoint);d.all_log_kp3b(n:endpoint)];
            end
            kpoly_individual_scatter = scatter(x,y,'filled',settings.shapes(LOOPY),'MarkerFaceColor', settings.colors(LOOPY),'MarkerEdgeColor','k');
            hold on
            n=n+PRMnum;
        end
    end
    
    xlabel(xlab)
    ylabel('log_2(kpoly)')
    
    legend(d.all_fh1_names_nobind, 'Location','northeastoutside');
    
     title(append(title_,' (',type,')'))
     subtitle(settings.settings_variable) 
     
    figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
    
    close all
    end

end