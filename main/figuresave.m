function figuresave(fig,file,figname)
% FIGURESAVE  Save figures to main pdf and as .fig for wholepackage.
    %
    %   FIGURESAVE(fig,file,figname) saves fig as a .fig called figname
    %           (and a png) and appends fig to the .pdf file.
    %
    %   Appends pdf of figure to the specified file in the current
    %   directory
    %
    % 
    %   The first time called, creates a new folder "figures" in the
    %   current directory and saves all future .fig and .png files there
    %   
    %   Inputs:
    %       fig     : figure to save
    %       file    : name of .pdf file to append figure to
    %       figname : .fig file name to save as 
    %
    %   See also MAKE_FORMIN_PLOTS, MAKE_FILAMENT_SCHEMATIC, WHOLEPACKAGE.
ax=fig;
persistent i; 
if isempty(i)
    saveas(ax,file);
    mkdir("figures");
    i = 0;
else
    exportgraphics(ax,file,'Append',true);
end
i = i+1;
savefig(ax,fullfile("figures",figname));
saveas(ax,fullfile("figures",append(extractBefore(figname,strlength(figname)-3),'.png')));
end
