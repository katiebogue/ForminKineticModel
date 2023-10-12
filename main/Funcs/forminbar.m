function fig=forminbar(forminList,settings,save)
% FORMINBAR  Makes and saves overview bargraphs for wholepackage.
    %
    %   Makes 3 bargraphs based on input structures (data and settings):
    %       1. Single, double, and dimer side-by-side for each formin
    %       2. Change upon NTD for each formin
    %       3. Number of PRMs for each formin (with PRMs)
    %   
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %   
    %   See also .

    arguments
        forminList Formin
        settings Options
        save logical=false
    end

    if save
        set(groot,'defaultfigureposition',[400 250 900 750])
    end

    fig(1)=figure;

    kpoly=[forminList.kpoly];
    fh1names=[forminList.name];

    % bar graph showing all polymerization rates
    kpoly_table = table([kpoly.single]', [kpoly.double]', [kpoly.dimer]', 'RowNames', fh1names);
   
    kpoly_bar = bar(kpoly_table{:,:});
    set(gca,'xtick',1:length(fh1names), 'xticklabel',kpoly_table.Properties.RowNames);
    xtickangle(90);
    maxBar = max(cellfun(@max, get(kpoly_bar, 'YData'))); 
    ylim([0, maxBar*1.1])
    set(kpoly_bar(1), 'FaceColor','b');
    set(kpoly_bar(2), 'FaceColor','r');
    set(kpoly_bar(3), 'FaceColor','g');
    legend( 'single', 'double', 'N-term dimerized');
    xlabel('Formins');
    ylabel('(kpoly)');
    
    title('Polymerization Rates per Formin')

    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end

    
    % Change in Polymerization Rates w/ Dimerization per Formin
    fig(2)=figure;
    kpoly_table_ratio = table([kpoly.ratio]', 'RowNames', fh1names);

    
    kpoly_bar_ratio = bar(kpoly_table_ratio{:,:});
    set(gca,'xtick',1:length(fh1names), 'xticklabel',kpoly_table_ratio.Properties.RowNames)
    xtickangle(90)
    maxBar = max(get(kpoly_bar_ratio, 'YData')); 
    ylim([0, maxBar*1.1])
    set(kpoly_bar_ratio(1), 'FaceColor','m')
    xlabel('Formins')
    ylabel('(kpoly N terminal dimerized/kpoly double)')

    title('Change in Polymerization Rates with NTD per Formin')
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end
    
    % Number of PRMs
    fig(3)=figure;
    hb = bar([forminList.PRMCount], 'stacked');
    set(hb(1), 'FaceColor','b');
    hold on
    
    xTick=get(gca,'xtick');
    set(gca,'xtick',1:length(fh1names), 'xticklabel', fh1names)
    xtickangle(90)
    maxBar = max(get(hb, 'YData')); 
    ylim([0, maxBar*1.1])
    
    xlabel('Formins')
    ylabel('Binding sites')
    title('Number of Binding Sites per Formin')
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end
end