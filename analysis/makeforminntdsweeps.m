function makeforminntdsweeps(lt,resultsloc)
% MAKEFORMINNTDSWEEPS generates Options, and Experiment for
% formins and saves heatmaps of simulated kpoly values accross multiple NTD
% locations
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

    %% Input files/ paths
    Nmax=400;
    
    pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
    forminfile="ForminTypes_knownG.txt"; % file containing sequences
    gatingfile="ForminTypes_gating_knownG.txt";

    %% create options object
    % modify this line to change the rate constants:
    opts=Options(lt,pythonpath,...
        "3st",...       % kpoly type
        resultsloc,...
        10^3.2778,...          % k_cap
        10^-2.8678,... % k_del
        10^3.15,...    % r_cap
        1,...    % r_del
        1);     % k_rel
    
    opts_1.r_cap_exp=0.57121;
    opts.set_equation(1); % using preset #1 (see Options class)
    opts.NTopt=2; 
    
    %% create experiment object
    Experiment1=Experiment(opts,forminfile,"uniprot",0.88); % using the input sequence option and concentration of profilin actin of 0.88
    Experiment1.set_gating_file(gatingfile);
    %Experiment1.set_gating("BNI1",0.5);
    
    %% Make heatmaps
    set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
    titles="log_2(k_{poly} N terminal dimerized/k_{poly} double)";
    
    NTDtable=makeNTDtable(Experiment1,20,Nmax);
    NTDtable2=makeNTDtable(Experiment1,1,Nmax);
    
    makefamilycurves(NTDtable)
    opts.update_results_folder
    opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_curves_","3st");
    figuresave(gcf,opts,append('NTDsweep_curves_','3st','.fig'),true);
    xlim([100 400])
    ylim([-0.2 0.2])
    figuresave(gcf,opts,append('NTDsweep_curves_truncated_','3st','.fig'),true);
    
    makefamilycurves(NTDtable2)
    figuresave(gcf,opts,append('NTDsweep_curves_noavg_','3st','.fig'),true);
    xlim([100 400])
    ylim([-0.3 0.3])
    figuresave(gcf,opts,append('NTDsweep_curves_truncated_noavg_','3st','.fig'),true);
    
    figure
    h1=makeheatmap(NTDtable,Nmax);
    h1.Title = {titles,"3st"};
    opts.update_results_folder
    opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","3st");
    figuresave(gcf,opts,append('NTDsweep_','3st','.fig'),true);
    
    figure
    h1=makeheatmap(NTDtable2,Nmax);
    h1.Title = {titles,"3st"};
    figuresave(gcf,opts,append('NTDsweep_noavg_','3st','.fig'),true);
    
    %%
    figure
    NTDmaxmin2=max_min_table(NTDtable2);
    b=make_barplot(NTDmaxmin2);
    title([titles,"3st"]);
    %figuresave(gcf,opts,append('NTDsweep_bar__noavg','3st','.fig'),true);
    %figuresave(gcf,opts,append('NTDsweep_bar__noavg','3st','.png'),true);
    exportgraphics(gcf,"NTDsweep_bar_3st.eps",'BackgroundColor','none','ContentType','vector')
    
    figure
    b=make_zerovis(NTDmaxmin2);
    title(["NT dist where kpoly ratio crosses 0","3st"]);
    %figuresave(gcf,opts,append('NTDsweep_zerobar__noavg','3st','.fig'),true);
    %figuresave(gcf,opts,append('NTDsweep_zerobar__noavg','3st','.png'),true);
    exportgraphics(gcf,"NTDsweep_zerobar_3st.eps",'BackgroundColor','none','ContentType','vector')
    
    %%
    figure
    NTDmaxmin=max_min_table(NTDtable);
    b=make_barplot(NTDmaxmin);
    title([titles,"3st"]);
    figuresave(gcf,opts,append('NTDsweep_bar__yesavg','3st','.fig'),true);
    
    figure
    b=make_zerovis(NTDmaxmin);
    title(["NT dist where kpoly ratio crosses 0","3st"]);
    figuresave(gcf,opts,append('NTDsweep_zerobar__yesavg','3st','.fig'),true);
end
function NTDtable=makeNTDtable(exp,rollingavg,Nmax)
% MAKENTDTABLE runs NTD_predictions for each formin and logs resulting
% kpoly values, ratios, and other stats in a table
arguments
    exp %Experiment with all of the formins to make predictions for
    rollingavg %number of nearby values to use in rolling avarge of ratios
    Nmax %max fh1 length value
end
    numformins=length(exp.ForminList);
    NTDtable=table('Size',[numformins*Nmax 10],'VariableTypes',["double","double","double","double","double","double","double","double","double","string"],'VariableNames',{'doubles','dimers','ratios','ratios_raw','NTD_dists','rollingstd','fh1length','std_dimer','std_double','formin_name'});
    for i=1:numformins
        formini=exp.ForminList(i);
        [doubles,dimers,ratios,NTD_dists,fh1ength]=NTD_predictions(formini,Nmax);
        x=(Nmax*(i-1)+1);
        y=Nmax*i;
        NTDtable(x:y,1)=array2table(doubles');
        NTDtable(x:y,2)=array2table(dimers');
        ratioavg=(movmean(ratios',rollingavg,"omitnan"));
        ratioavg(isnan(ratios'))=NaN;
        NTDtable(x:y,3)=array2table(ratioavg);
        NTDtable(x:y,4)=array2table(ratios');
        NTDtable(x:y,5)=array2table(NTD_dists');
        NTDtable(x:y,6)=array2table(movstd(ratios',5,"omitnan"));
        NTDtable(x:y,7)=array2table(fh1ength');
        NTDtable(x:y,8)=array2table((movstd(dimers',5,"omitnan"))./mean(dimers,"omitnan"));
        NTDtable(x:y,9)=array2table((movstd(doubles',5,"omitnan"))./mean(doubles,"omitnan"));
        NTDtable(x:y,10)={formini.name};
    end
end
function [doubles,dimers,ratios,NTD_dists,fh1ength]=NTD_predictions(formin,Nmax)
% NTD_predictions sweeps through all possible NTD locations up to fh1 length=Nmax for
% the input formin and returns kpoly values, corresponding NTD distances
% and FH1 lengths
arguments
    formin %formin to make predicitons for
    Nmax %max fh1 length value
end
    doubles=NaN(1,Nmax);
    dimers=NaN(1,Nmax);
    ratios=NaN(1,Nmax);
    NTD_dists=[1:Nmax];
    fh1ength=NaN(1,Nmax);
    for i=0:(Nmax-1)
        if formin.length<Nmax+1
            kpoly=formin.kpoly;
            NTD_dist=formin.PRMList(1,formin.PRMCount).dist_NT;
            doubles(NTD_dist)=kpoly.double;
            dimers(NTD_dist)=kpoly.dimer;
            ratioval=log2(kpoly.ratio);
            fh1ength(NTD_dist)=formin.length;
            if ratioval==-Inf
                ratioval=NaN;
            end
            ratios(NTD_dist)=ratioval;
            formin.add_length(1)
        else
            formin.add_length(-i)
            return
        end
    end
end

function makefamilycurves(tab)
    formins=unique(tab.formin_name);
    colors=makepoints(length(formins)+2,'w');
    [B,I]=sort(arrayfun(@(x) hex2dec(erase(x(1), "#")), colors));
    colors=colors(I);
    linestyles=["-","--","-."];
    linetype=1;
    figure
    hold on
    for i=1:length(formins)
        forminname=formins(i);
        subtab= tab(strcmp(tab.formin_name, forminname), :);
        plot(subtab.NTD_dists,subtab.ratios,"LineWidth",2.5,Color=colors(i+1),LineStyle=linestyles(linetype))
        if linetype==3
            linetype=1;
        else
            linetype=linetype+1;
        end
    end

    xlabel('NTD dist')
    ylabel('log_2 Polymerization Rate Ratio (dimer/double)')

    legend(formins)
    hold off
end
function h=makeheatmap(tab,Nmax)
% MAKEHEATMAP creates heatmap of Kpoly ratios for each formin for each
% possible NTD dist, from a table generated by makeNTDtable
arguments 
    tab %table generated by makeNTDtable
    Nmax %max fh1 length value
end
    %ogtab=tab;
    % x=movmean(tab.ratios,10,"omitnan");
    % tab.ratios=x;
    h = heatmap(tab,'formin_name','NTD_dists','ColorVariable','ratios');
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal';
    load('customcolorbar_red_blue_large.mat');
    h.Colormap=CustomColormap;
    max_min_tab=max_min_table(tab);
    allmaxmin=[max_min_tab.max; abs(max_min_tab.min)];
    allmaxmin=sort(allmaxmin);
    allmaxmin = allmaxmin(~isnan(allmaxmin));
    maxratio=allmaxmin(end);
    if maxratio-allmaxmin(end-1)>5
        maxratio=allmaxmin(end-1);
    end
    h.ColorLimits=[-maxratio,maxratio];
    h.ColorLimits=[-5,5];
    yvals=[1:Nmax];
    CustomYLabels = string(yvals);
    CustomYLabels(mod(yvals,20) ~= 0) = " ";
    h.YDisplayLabels = CustomYLabels;
    h.GridVisible = 'off';

end
function h=makeheatmap_std(tab,Nmax,disttype,stdtype)
% MAKEHEATMAP_STD creates heatmap of standard deviations of kpoly/ratios 
% for each formin for each possible NTD dist, from a table generated by makeNTDtable
arguments 
    tab %table generated by makeNTDtable
    Nmax %max fh1 length value
    disttype % y axis for heatmaps, either NTD distances or total FH1 length
    stdtype % color variable for heatmaps, a column in tab that contains std values
end
    h = heatmap(tab,'formin_name',disttype,'ColorVariable',stdtype);
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal';
    tab.ratios=tab.(stdtype);
    max_min_tab=max_min_table(tab);
    allmaxmin=[max_min_tab.max; abs(max_min_tab.min)];
    allmaxmin=sort(allmaxmin);
    allmaxmin = allmaxmin(~isnan(allmaxmin));
    maxratio=allmaxmin(end);
    if maxratio-allmaxmin(end-1)>1
        maxratio=allmaxmin(end-1);
    end
    h.ColorLimits=[0,maxratio];
    yvals=str2double(h.YData);
    CustomYLabels = string(yvals);
    CustomYLabels(mod(yvals,20) ~= 0) = " ";
    h.YDisplayLabels = CustomYLabels;
    h.GridVisible = 'off';
end
function mmtab=max_min_table(tab)
% MAX_MIN_TABLE creates table of max and min ratios for each formin in the
% input table
    formins=unique(tab.formin_name);
    numformins=length(formins);
    mmtab=table('Size',[numformins 4],'VariableTypes',["double","double","string","double"],'VariableNames',{'max','min','formin_name','zeropt'});
    mmtab.formin_name=formins;
    for i=1:numformins
        formin_name=mmtab.formin_name(i);
        [ext,zeropt]=get_max_min(tab,formin_name);
        mmtab{i,'max'}=ext(1);
        mmtab{i,'min'}=ext(2);
        mmtab{i,'zeropt'}=zeropt;
    end
end
function [ext,zeropt]=get_max_min(tab,formin_name)
% GET_MAX_MIN finds max and min ratios for a given formin found in the
% input table
    rows=matches(tab.formin_name,formin_name);
    vals=tab{rows,'ratios'};
    zeropt= find(vals >= 0,1);
    if isempty(zeropt)
        zeropt=0;
    end
    ext=[max(vals), min(vals)];
end 
function b=make_barplot(mmtab)
% MAKE_BARPLOT creates a barplot of the max and min values in the input
% mmtab, which should be generated by max_min_tab
    b=bar(categorical(mmtab.('formin_name')'),[mmtab.('min')';mmtab.('max')'],'stacked');
    xtickangle(90)
    legend({'minimum','maximum'})
    b(1).FaceColor='blue';
    b(2).FaceColor='red';
end

function b=make_zerovis(mmtab)
% 
    b=bar(categorical(mmtab.('formin_name')'),[mmtab.('zeropt')']);
    xtickangle(90)
    b(1).FaceColor='blue';
end