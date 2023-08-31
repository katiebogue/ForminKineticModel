%makes all of the correlation plot figures at the end of wholepackage

%% Overview figures
make_bargraph_overviews(Data,settings)

%% Formin correlation plots

% Change in Polymerization Rates vs Number of PRMs
NTDplot(all_iSite_tot,'Number of PRMs','Change in Polymerization Rates vs Number of PRMs',Data,settings);

% Polymerization Rates vs Number of PRMs
kpolyplot(all_iSite_tot,'Number of PRMs',20,Data,settings);

% Change in Polymerization Rates vs Length of FH1 Domain
NTDplot(all_fh1_length,'Length of FH1 domain (1st PRM to FH2)','Change in Polymerization Rates vs Length of FH1 Domain',Data,settings);

% Polymerization Rates vs Length of FH1 Domain
kpolyplot(all_fh1_length,'Length of FH1 domain (1st PRM to FH2)',400,Data,settings);

% Change in Polymerization Rates vs Mean PRM size
NTDplot(all_mean_PP_length,'Mean PRM size','Change in Polymerization Rates vs Mean PRM size',Data, settings);

% Change in Polymerization Rates vs Mean PRM size x Number of PRMs
NTDplot(all_PP_length_x_PP_isite_tot,'Mean PRM size x #PRM','Change in Polymerization Rates vs Mean PRM size x Number of PRMs',Data, settings);

%% Per PRM plots

% Get colors and shapes for plots
[settings.colors,settings.shapes]=makepoints(length(all_fh1_names_nobind),'w');

% Vs. PRM length
PRMplot(all_PP_length,'Length of PRM',Data,settings);

% Distance from PRM to FH2
PRMplot(all_PP_location,'Distance from PRM to FH2',Data,settings);

% Distance from PRM to end
PRMplot(all_PP_dist_end,'Distance from PRM to N-term',Data,settings);

% Change in Polymerization Rates vs. PP length per individual PRM
PRM_NTD_plot(all_PP_length,'Length of PRM',Data,settings)

% Change in Polymerization Rates vs. PP dist to FH2 per individual PRM
PRM_NTD_plot(all_PP_location,'Distance from PRM to FH2',Data,settings)

% Change in Polymerization Rates vs. PP dist to NT per individual PRM
PRM_NTD_plot(all_PP_dist_end,'Distance from PRM to N-term',Data,settings)

%% Polymer stats plots

% Change in polymer stats vs. FH1 length per individual PRM
polymerstat_change_plot(all_fh1_length_PP,'FH1 length',Data,settings);

% Change in polymer stats vs. PP dist to NT per individual PRM
polymerstat_change_plot(all_PP_dist_end,'Distance from PRM to N-term',Data,settings);

% Change in polymer stats vs. Fractional distance from FH2 per individual PRM
all_frac_dist=all_PP_location./all_fh1_length_PP;
polymerstat_change_plot(all_frac_dist,'Fractional Distance from PRM to FH2',Data,settings);

%% Make overview tables
make_formin_tables



