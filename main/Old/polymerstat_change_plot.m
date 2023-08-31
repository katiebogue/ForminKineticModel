function polymerstat_change_plot(parameter,xlab,d,settings)
% POLYMERSTAT_CHANGE_PLOT  Creates and saves a scatterplot of change in polymer
%                  statistics vs. a parameter.
    %
    %   POLYMERSTAT_CHANGE_PLOT(x,'xlab',Data,settings) creates and saves a
    %   scatterplot of change in polymer statistics vs. x. 
    %   
    %   Plots change in 3 polymer stats:
    %       1. pr_0
    %       2. 1-pocc
    %       3. 1-pocc_0
    %
    %   Inputs:
    %       parameter   : array containing values to plot against change in
    %                     polymer statistics
    %       xlab        : x-axis (parameter) label (string)
    %       d           : structure containing data
    %                     must include: 
    %                     all_p_occ2a, all_p_occ2a_0, all_p_r2a,
    %                     all_p_occ2b, all_p_occ2b_0, all_p_r2b,
    %                     all_p_occ3a, all_p_occ3a_0, all_p_r3a,
    %                     all_p_occ3b, all_p_occ3b_0, all_p_r3b, 
    %                     all_iSite_tot, all_fh1_names_nobind
    %       settings    : structure containing settings variables
    %                     must include: 
    %                     settings_variable, pdf_name, shapes, colors
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS, FIND_POCC, WHOLEPACKAGE.

    ya=d.all_p_r3a ./ d.all_p_r2a;
    yb=d.all_p_r3b ./ d.all_p_r2b;
    plotchange(ya,yb,'Concentration at the Barbed End (pr_{0})');

    ya=(1-d.all_p_occ3a) ./ (1-d.all_p_occ2a);
    yb=(1-d.all_p_occ3b) ./ (1-d.all_p_occ2b);
    plotchange(ya,yb,'PRM Non-occlusion Probability (1-pocc)');

    ya=(1-d.all_p_occ3a_0) ./ (1-d.all_p_occ2a_0);
    yb=(1-d.all_p_occ3b_0) ./ (1-d.all_p_occ2b_0);
    plotchange(ya,yb,'Barbed End Non-occlusion Probability (1-pocc_{0})');
    
    function plotchange(ya,yb,ytxt)
        % ya is 3a/2a, yb is 3b/2b
        ylab=append(ytxt,' Dimer/Double');
        title_=({append('Change in ',ytxt),append(' with NTD vs. ',xlab,' per individual PRM')});
        titlefig=(append('Change in ',ytxt,' with NTD vs. ',xlab,' per individual PRM'));

        %% Plot based on each filament
        h(1)=figure;
        scatter_3a_2a = scatter(parameter,ya,'filled', 'o');
        hold on
        scatter_3b_2b = scatter(parameter,yb, 'filled', 's');
        hold on
        
        xlabel(xlab)
        ylabel(ylab)
        legend('filament a', 'filament b');
        
        title(title_)
        subtitle(settings.settings_variable)
         
        figuresave(gcf,settings.pdf_name,append('stats.',titlefig,'.fig'));
        
        
        %% Plot based on formin
        h(2)=figure;
        n=1;
        
        while n<= length(d.all_p_r2a)
            for LOOPY = 1:length(d.all_iSite_tot)
                PRMnum = d.all_iSite_tot(LOOPY);
                endpoint = n + PRMnum-1;
                x=[parameter(n:endpoint);parameter(n:endpoint)];
                y=[ya(n:endpoint);yb(n:endpoint)];
                scatter_formin = scatter(x,y, 'filled',settings.shapes(LOOPY), 'MarkerFaceColor', settings.colors(LOOPY),'MarkerEdgeColor','k');
                hold on
                n=n+PRMnum;
            end
        end
        
        xlabel(xlab)
        ylabel(ylab)
        
        legend(d.all_fh1_names_nobind, 'Location','southoutside','NumColumns',6,'FontSize',8);
        
        title(title_)
        subtitle(settings.settings_variable) 
         
        figuresave(gcf,settings.pdf_name,append('stats.',titlefig,'.fig'));
        savefig(h,fullfile("figures",append('stats.',titlefig,'.fig')));
        
        close all
    end
end