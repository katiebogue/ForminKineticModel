% MAKE_FORMIN_TABLES Makes and saves uitable .figs for per formin and per
%                    PRM data for wholepackage.
    %
    %   Makes 2 uitables and saves them as .figs:
    %       1. Data on a per formin basis
    %       2. Data on a per PRM basis
    %
    %   Different variables are included based on variable delivery
    %   
    %   See also MAKE_FORMIN_PLOTS, WHOLEPACKAGE. 
    
%% Per formin data

fig= figure('Name','Per Formin Data');
if delivery == 'Y'
    uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto',...
        'ColumnName',{'k_poly Single (log2)', 'k_poly Double (log2)', 'k_poly N-Dimer (log2)','log_2(kpoly ratio)','# Binding Sites','FH1 length', 'Mean PRM size','Mean PRM size x # Binding Sites','Single Capture Rate (log2)','Double Capture Rate (log2)','Dimer Capture Rate (log2)','Single Delivery Rate (log2)','Double Delivery Rate (log2)','Dimer Delivery Rate (log2)','Profilin-Actin Concentration|μM|c_PA','Capture rate constant|μM^(-1)s^(-1)|k_paf','Delivery rate constant|μM^(-1)s^(-1)|k_pab','Reverse Capture Rate|s^(-1)|k_paf_rev','Ring opening Rate|s^(-1)|r_PF_rev','Reverse Delivery Rate|s^(-1)|r_paf_rev'},...
        'Data',[all_kpoly1_nobind, all_kpoly2_nobind, all_kpoly3_nobind, all_log_kpoly3_2_nobind, all_iSite_tot, all_fh1_length, all_mean_PP_length, all_PP_length_x_PP_isite_tot, all_kcap1_nobind, all_kcap2_nobind, all_kcap3_nobind,all_kdel1_nobind,all_kdel2_nobind,all_kdel3_nobind,ones(length(all_fh1_names_nobind),1).*c_PA,ones(length(all_fh1_names_nobind),1).*k_paf,ones(length(all_fh1_names_nobind),1).*k_pab,ones(length(all_fh1_names_nobind),1).*k_paf_rev,ones(length(all_fh1_names_nobind),1).*r_PF_rev,ones(length(all_fh1_names_nobind),1).*r_paf_rev]);
end
if delivery == 'N'
    uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto',...
        'ColumnName',{'k_poly Single (log2)', 'k_poly Double (log2)', 'k_poly N-Dimer (log2)','log_2(kpoly ratio)','# Binding Sites','FH1 length', 'Mean PRM size','Mean PRM size x # Binding Sites','Profilin-Actin Concentration|μM|c_PA','Capture rate constant|μM^(-1)s^(-1)|k_paf','Delivery rate constant|μM^(-1)s^(-1)|k_pab','Reverse Capture Rate|s^(-1)|k_paf_rev'},...
        'Data',[all_kpoly1_nobind, all_kpoly2_nobind, all_kpoly3_nobind, all_log_kpoly3_2_nobind, all_iSite_tot, all_fh1_length, all_mean_PP_length, all_PP_length_x_PP_isite_tot, ones(length(all_fh1_names_nobind),1).*c_PA,ones(length(all_fh1_names_nobind),1).*k_paf,ones(length(all_fh1_names_nobind),1).*k_pab,ones(length(all_fh1_names_nobind),1).*k_paf_rev]);
end

uit.RowName = all_fh1_names_nobind;

saveas(fig,'Per Formin Data_'+fig_name);


close all
%% Per PRM data 
fig= figure('Name','Per PRM Data');
if delivery == 'Y'
    uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto',...
        'ColumnName',{'distance from PRM to FH2','PRM length','FH1 Length','Single Capture','Double Capture a','Double Capture b','Dimer Capture a','Dimer Capture b','log_2(Capture Dimer/Double) a','log_2(Capture Dimer/Double) b','Single Delivery','Double Delivery a','Double Delivery b','Dimer Delivery a','Dimer Delivery b','log_2((Delivery Dimer/Double) a','log_2(Delivery Dimer/Double) b','Single KPoly','Double KPoly a','Double KPoly b','Dimer KPoly a','Dimer KPoly b','log_2(KPoly Dimer/Double) a','log_2(KPoly Dimer/Double) b','Profilin-Actin Concentration|μM|c_PA','Capture rate constant|μM^(-1)s^(-1)|k_paf','Delivery rate constant|μM^(-1)s^(-1)|k_pab','Reverse Capture Rate|s^(-1)|k_paf_rev','Ring opening Rate|s^(-1)|r_PF_rev','Reverse Delivery Rate|s^(-1)|r_paf_rev','Single Pocc|p_occ1','Double Pocc a|p_occ2a','Double Pocc b|p_occ2b','Dimer Pocc a|p_occ3a','Dimer Pocc b|p_occ3b','Single Pocc_0|p_occ1_0','Double Pocc_0 a|p_occ2a_0','Double Pocc_0 b|p_occ2b_0','Dimer Pocc_0 a|p_occ3a_0','Dimer Pocc_0 b|p_occ3b_0','Single PRM barbed end concentration|p_r1','Double PRM barbed end concentration a|p_r2a','Double PRM barbed end concentration b|p_r2b','Dimer PRM barbed end concentration a|p_r3a','Dimer PRM barbed end concentration b|p_r3b'},...
        'Data',[all_PP_location,all_PP_length,all_fh1_length_PP,all_kc1,all_kc2a,all_kc2b,all_kc3a,all_kc3b,all_kcap3a_2a,all_kcap3b_2b,all_kd1,all_kd2a,all_kd2b,all_kd3a,all_kd3b,all_kdel3a_2a,all_kdel3b_2b,all_kp1,all_kp2a,all_kp2b,all_kp3a,all_kp3b,all_kpoly3a_2a,all_kpoly3b_2b,ones(length(all_fh1_names_perPRM),1).*c_PA,ones(length(all_fh1_names_perPRM),1).*k_paf,ones(length(all_fh1_names_perPRM),1).*k_pab,ones(length(all_fh1_names_perPRM),1).*k_paf_rev,ones(length(all_fh1_names_perPRM),1).*r_PF_rev,ones(length(all_fh1_names_perPRM),1).*r_paf_rev,all_p_occ1,all_p_occ2a,all_p_occ2b,all_p_occ3a,all_p_occ3b,all_p_occ1_0,all_p_occ2a_0,all_p_occ2b_0,all_p_occ3a_0,all_p_occ3b_0,all_p_r1,all_p_r2a,all_p_r2b,all_p_r3a,all_p_r3b]);
end
if delivery == 'N'
    uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto',...
        'ColumnName',{'distance from PRM to FH2','PRM length','FH1 Length','Single KPoly','Double KPoly a','Double KPoly b','Dimer KPoly a','Dimer KPoly b','log_2(KPoly Dimer/Double) a','log_2(KPoly Dimer/Double) b','Profilin-Actin Concentration|μM|c_PA','Capture rate constant|μM^(-1)s^(-1)|k_paf','Delivery rate constant|μM^(-1)s^(-1)|k_pab','Reverse Capture Rate|s^(-1)|k_paf_rev','Single Pocc|p_occ1','Double Pocc a|p_occ2a','Double Pocc b|p_occ2b','Dimer Pocc a|p_occ3a','Dimer Pocc b|p_occ3b','Single Pocc_0|p_occ1_0','Double Pocc_0 a|p_occ2a_0','Double Pocc_0 b|p_occ2b_0','Dimer Pocc_0 a|p_occ3a_0','Dimer Pocc_0 b|p_occ3b_0','Single PRM barbed end concentration|p_r1','Double PRM barbed end concentration a|p_r2a','Double PRM barbed end concentration b|p_r2b','Dimer PRM barbed end concentration a|p_r3a','Dimer PRM barbed end concentration b|p_r3b'},...
        'Data',[all_PP_location,all_PP_length,all_fh1_length_PP,all_kp1,all_kp2a,all_kp2b,all_kp3a,all_kp3b,all_kpoly3a_2a,all_kpoly3b_2b,ones(length(all_fh1_names_perPRM),1).*c_PA,ones(length(all_fh1_names_perPRM),1).*k_paf,ones(length(all_fh1_names_perPRM),1).*k_pab,ones(length(all_fh1_names_perPRM),1).*k_paf_rev,all_p_occ1,all_p_occ2a,all_p_occ2b,all_p_occ3a,all_p_occ3b,all_p_occ1_0,all_p_occ2a_0,all_p_occ2b_0,all_p_occ3a_0,all_p_occ3b_0,all_p_r1,all_p_r2a,all_p_r2b,all_p_r3a,all_p_r3b]);
end
uit.RowName = all_fh1_names_perPRM;

saveas(fig,'Per PRM Data_'+fig_name);

close all