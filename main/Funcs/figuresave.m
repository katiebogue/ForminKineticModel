function figuresave(fig,options,figname,notestab)
% FIGURESAVE  Save figures to main pdf and as .fig.
    %
    %   FIGURESAVE(fig,options,figname) saves fig as a .fig called figname
    %           (and a png) and appends fig to the .pdf file.
    %
    %   Appends pdf of figure to the file in the directory specified in
    %   options
    %
    % 
    %   Saves/appends figure to pdf "RESULTS_ date/time" in
    %   options.resultsdir/options.resultsfolder (or creates folder if
    %   nonexistant).
    %   Saves .fig and .png of figure to 
    %   options.resultsdir/options.resultsfolder/figures 
    %   (or creates folder if nonexistant).
    %   
    %   Inputs:
    %       fig     : figure to save
    %       options : Options class
    %       figname : .fig file name to save as 
    %       notestab: table to add to fist page if nothing has been saved
    %                 yet
    %
    %   See also .
   
    arguments
        fig
        options Options
        figname string
        notestab=-1
    end
    ax=fig;
    if ~exist(fullfile(options.resultsdir,options.resultsfolder),'dir')
        mkdir (options.resultsdir,options.resultsfolder)
        savesnotes();
        saveas(ax,fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"));
        mkdir(fullfile(options.resultsdir,options.resultsfolder),"figures");
    elseif ~exist(fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'file')
        savesnotes();
        saveas(ax,fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"));
        if ~exist(fullfile(options.resultsdir,options.resultsfolder,"figures"),'dir')
            mkdir(fullfile(options.resultsdir,options.resultsfolder),"figures");
        end
    else
        if ~exist(fullfile(options.resultsdir,options.resultsfolder,"figures"),'dir')
            mkdir(fullfile(options.resultsdir,options.resultsfolder),"figures");
        end
        exportgraphics(ax,fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'Append',true);
    end
    savefig(ax,fullfile(options.resultsdir,options.resultsfolder,"figures",figname));
    saveas(ax,fullfile(options.resultsdir,options.resultsfolder,"figures",append(extractBefore(figname,strlength(figname)-3),'.png')));

    function savesnotes()
        if notestab==-1
            return
        end
        import mlreportgen.dom.*;
        import mlreportgen.report.*;
        rpt=Report(fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'pdf');
        table=Table(notestab);
        slicer=mlreportgen.utils.TableSlicer("Table",table);
        slices=slicer.slice();
        ch=Chapter();
        add(ch,slices.Table);
        add(rpt,ch);
        close(rpt)
    end
end
