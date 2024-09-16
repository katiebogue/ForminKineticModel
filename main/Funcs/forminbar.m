function fig=forminbar(forminList,settings,save,NameValueArgs)
% FORMINBAR  Makes and saves overview bargraphs for a collection of
% formins.
    %   fig = FORMINBAR(forminList,settings) creates 3 bargraphs of per 
    %   formin information 
    %
    %   fig = FORMINBAR(forminList,settings,true) creates and saves 3 
    %   bargraphs of per formin information
    %
    %   fig = FORMINBAR(forminList,settings,true/false,kpolyscale='scale',ratioscale='scale') 
    %   creates and (if save=true) saves 3 bargraphs of per formin 
    %   information using the specified scaling
    %
    %
    %   Makes 3 bargraphs:
    %       1. Single, double, and dimer polymerization rates side-by-side 
    %          for each formin
    %       2. Change upon NTD for each formin
    %       3. Number of PRMs for each formin (with PRMs)
    %   
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %       kpolyscale  : y-axis scale (can be none, log2, log10, ln) for
    %                     kpoly bargrph (Name-value argument)
    %       ratioscale  : y-axis scale (can be none, log2, log10, ln) for
    %                     NTD change bargrph (Name-value argument)
    %
    %   Calls figuresave
    %   
    %   See also FIGURESAVE, FORMIN, PRM, EXPERIMENT, OPTIONS.

    arguments
        forminList Formin
        settings Options
        save logical=false
        NameValueArgs.kpolyscale string {mustBeMember(NameValueArgs.kpolyscale,{'none','log2','log10','ln'})}
        NameValueArgs.ratioscale string {mustBeMember(NameValueArgs.ratioscale,{'none','log2','log10','ln'})}
    end

  
    if isfield(NameValueArgs,"kpolyscale")
        scale=NameValueArgs.kpolyscale;
    else
        scale="none";
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

   
    kpoly=[forminList.kpoly];
    fh1names=[forminList.name];

    if length(fh1names)>1
        fig(1)=figure;
    
        % bar graph showing all polymerization rates
        kpoly_table = table(yscale([kpoly.single]'), yscale([kpoly.double]'), yscale([kpoly.dimer]'), 'RowNames', fh1names);
       
        kpoly_bar = bar(kpoly_table{:,:});
        set(gca,'xtick',1:length(fh1names), 'xticklabel',kpoly_table.Properties.RowNames);
        xtickangle(90);
        vals=get(kpoly_bar, 'YData');
        maxBar = max(cellfun(@max, vals)); 
        ylim([0, maxBar*1.1])
        set(kpoly_bar(1), 'FaceColor','b');
        set(kpoly_bar(2), 'FaceColor','r');
        set(kpoly_bar(3), 'FaceColor','g');
        legend( 'single', 'double', 'N-term dimerized');
        xlabel('Formins');
        ylabel(ylab);
        
        title('Polymerization Rates per Formin')
    
        if save
            figuresave(gcf,settings,append(gca().Title.String,'.fig'));
        end

    end
    
    % Change in Polymerization Rates w/ Dimerization per Formin

    if isfield(NameValueArgs,"ratioscale")
        scale=NameValueArgs.ratioscale;
    else
        scale="log2";
    end

    if scale=="none"
        yscale=@(x) x;
        ylab="k_{poly} N terminal dimerized/k_{poly} double";
    elseif scale=="log2"
        ylab="log_{2}(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log(x);
    end

    fig(2)=figure;
    kpoly_table_ratio = table(yscale([kpoly.ratio]'), 'RowNames', fh1names);

    
    kpoly_bar_ratio = bar(kpoly_table_ratio{:,:});
    set(gca,'xtick',1:length(fh1names), 'xticklabel',kpoly_table_ratio.Properties.RowNames)
    xtickangle(90)
    maxBar = max(get(kpoly_bar_ratio, 'YData')); 
    minBar=min([get(kpoly_bar_ratio, 'YData'),0]); 
    ylim([minBar*1.1, maxBar*1.1])
    set(kpoly_bar_ratio(1), 'FaceColor','m')
    xlabel('Formins')
    ylabel(ylab)

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