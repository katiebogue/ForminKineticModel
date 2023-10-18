function figuresave(fig,options,figname)
% FIGURESAVE  Save figures to main pdf and as individual .pngs and .figs.
    %
    %   FIGURESAVE(fig,options,figname) saves fig as a .fig and .png called 
    %   figname and appends fig to the .pdf file specified in options.
    %
    %   Does the following:
    %   - Saves/appends figure to pdf "RESULTS.pdf" in
    %     options.resultsdir/options.resultsfolder (or creates folder if
    %     nonexistant).
    %   - Saves .fig and .png of figure to 
    %     options.resultsdir/options.resultsfolder/figures 
    %     (or creates folder if nonexistant).
    %   - If RESULTS.pdf does not already exist, the pdf is created and the
    %     table in options.optionstable is printed to the first page of 
    %     the pdf.
    %   
    %   Inputs:
    %       fig     : figure to save
    %       options : Options class
    %       figname : .fig file name to save as 
    %
    %   See also OPTIONS.
   
    arguments
        fig
        options Options
        figname string
    end
    ax=fig;
    if ~exist(fullfile(options.resultsdir,options.resultsfolder),'dir')
        mkdir (options.resultsdir,options.resultsfolder)
        savesnotes();
        exportgraphics(ax,fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'Append',true);
        mkdir(fullfile(options.resultsdir,options.resultsfolder),"figures");
    elseif ~exist(fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'file')
        savesnotes();
        exportgraphics(ax,fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'Append',true);
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
        import mlreportgen.dom.*;
        import mlreportgen.report.*;
        rpt=Report(fullfile(options.resultsdir,options.resultsfolder,"RESULTS.pdf"),'pdf');
        table=Table(options.optionstable);
        slicer=mlreportgen.utils.TableSlicer("Table",table);
        slices=slicer.slice();
        ch=Chapter();
        add(ch,slices.Table);
        add(rpt,ch);
        close(rpt)
    end
end
