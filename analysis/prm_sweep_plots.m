% PRM_SWEEP_PLOTS makes various scatterplots from PRM sweep table
% 
% Loads sweep_workspace.mat and makes the following scatterplots (using
% gscatter, so grouped as specified) from the variable PRM_sweep_table:
%   1- Max Kpoly ratio vs. PRM location; grouped by FH1 length
%   2- Min Kpoly ratio vs. PRM location; grouped by FH1 length 
%   3- PRM size where the max Kpoly ratio is vs. PRM location; grouped by FH1 length 
%   4- PRM size where the min Kpoly ratio is vs. PRM location; grouped by FH1 length 
%   5- PRM size of boundry of max Kpoly ratio (within 0.05xstd of max) vs. PRM location; grouped by FH1 length  
%   6- PRM size of boundry of min Kpoly ratio (within 0.05xstd of max) vs. PRM location; grouped by FH1 length   
%   7- PRM size of max-min Kpoly ratio vs. PRM location; grouped by FH1 length   
%   8- PRM size of max-min Kpoly ratio vs. PRM location; grouped by model   
%   9- PRM size of max-min Kpoly ratio vs. FH1 length; grouped by PRM location  
% 
% See also DISTINGUISHABLE_COLORS, GSCATTER, FIGURESAVE, OPTIONS,
% GENERATE_SWEEP_DATA, PRM_SIZE_SWEEP_INFO

load("sweep_workspace.mat")

colors=distinguishable_colors(40,'w');

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.max_val,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','Max Kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("Max Kpoly ratios for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.min_val,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','Min Kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("Min Kpoly ratios for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.max_x,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','PRM size of max kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("PRM size of max kpoly ratio for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.min_x,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','PRM size of min kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("PRM size of min kpoly ratio for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.max_x_boundry,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','PRM size of boundry of max kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("PRM size of boundry of max kpoly ratio (within 0.05xstd of max) for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.min_x_boundry,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','PRM size of boundry of min kpoly ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("PRM size of boundry of min kpoly ratio (within 0.05xstd of min) for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.max_x-PRM_sweep_table.min_x,PRM_sweep_table.FH1_length,colors,'.',8,'on','Distance from PRM to FH2','PRM size of max ratio - min ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'FH1 length';
title1=strcat("PRM size of max ratio - min ratio for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.PRM_loc,PRM_sweep_table.max_x-PRM_sweep_table.min_x,PRM_sweep_table.model,colors,'.',8,'on','Distance from PRM to FH2','PRM size of max ratio - min ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'model';
title1=strcat("PRM size of max ratio - min ratio for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))

figure
h=gscatter(PRM_sweep_table.FH1_length,PRM_sweep_table.max_x-PRM_sweep_table.min_x,PRM_sweep_table.PRM_loc,colors,'.',8,'on','FH1 length','PRM size of max ratio - min ratio');
lgd=legend('Location','northeastoutside');
lgd.Title.String = 'Distance from PRM to FH2';
title1=strcat("PRM size of max ratio - min ratio for PRM size sweep");
title(title1)
figuresave(gcf,Opts2,strcat(title1,".fig"))
