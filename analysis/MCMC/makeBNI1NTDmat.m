function [bni1mat, h, NTD_dists,fh1ength] = makeBNI1NTDmat(mcmcfile)

    mcmcresultsloc=fullfile(mcmcfile,"mcmc_results.mat");

    % if isempty(who('-file',mcmcfile,'vals_trueparamconvert'))
    %     vals_trueparamconvert = trueparamconvert(mcmcfile);
    %     save(mcmcfile, 'vals_trueparamconvert',"-append");
    % else
    %     load(mcmcfile, 'vals_trueparamconvert');
    % end
    % 
    % 
    % params=10.^(vals_trueparamconvert);
    % params(4)=vals_trueparamconvert(4); %rcap_exp
    % 
    % clear vals_trueparamconvert

    load("forminexperimentobjs.mat","Experiment_BNI1")
    curlength=Experiment_BNI1.ForminList.length;
    lenchange=106-curlength;
    Experiment_BNI1.ForminList.add_length(lenchange);
    % formin1=Experiment_BNI1.ForminList;
    % Nmax=450;
    % Nmin=106;

    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, mcmcresultsloc,0, 0, 0);
    NTDdist=Experiment_BNI1.ForminList.PRMList(1,Experiment_BNI1.ForminList.PRMCount).dist_NT;
    sgtitle(strcat("BNI1; FH1 length: ", num2str(Experiment_BNI1.ForminList.length)," NT dist: ", num2str(NTDdist)))
    %saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'epsc');
    exportgraphics(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'BackgroundColor','none','ContentType','vector');
    saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.png')),'png');

    Experiment_BNI1.ForminList.add_length(49);
    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, mcmcresultsloc,0, 0, 0);
    NTDdist=Experiment_BNI1.ForminList.PRMList(1,Experiment_BNI1.ForminList.PRMCount).dist_NT;
    sgtitle(strcat("BNI1; FH1 length: ", num2str(Experiment_BNI1.ForminList.length)," NT dist: ", num2str(NTDdist)))
    %saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'epsc');
    exportgraphics(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'BackgroundColor','none','ContentType','vector');
    saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.png')),'png');

    Experiment_BNI1.ForminList.add_length(50);
    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, mcmcresultsloc,0, 0, 0);
    NTDdist=Experiment_BNI1.ForminList.PRMList(1,Experiment_BNI1.ForminList.PRMCount).dist_NT;
    sgtitle(strcat("BNI1; FH1 length: ", num2str(Experiment_BNI1.ForminList.length)," NT dist: ", num2str(NTDdist)))
    %saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'epsc');
    exportgraphics(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'BackgroundColor','none','ContentType','vector');
    saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.png')),'png');

    Experiment_BNI1.ForminList.add_length(100);
    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, mcmcresultsloc,0, 0, 0);
    NTDdist=Experiment_BNI1.ForminList.PRMList(1,Experiment_BNI1.ForminList.PRMCount).dist_NT;
    sgtitle(strcat("BNI1; FH1 length: ", num2str(Experiment_BNI1.ForminList.length)," NT dist: ", num2str(NTDdist)))
    %saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'epsc');
    exportgraphics(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'BackgroundColor','none','ContentType','vector');
    saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.png')),'png');

    Experiment_BNI1.ForminList.add_length(200);
    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, mcmcresultsloc,0, 0, 0);
    NTDdist=Experiment_BNI1.ForminList.PRMList(1,Experiment_BNI1.ForminList.PRMCount).dist_NT;
    sgtitle(strcat("BNI1; FH1 length: ", num2str(Experiment_BNI1.ForminList.length)," NT dist: ", num2str(NTDdist)))
    %saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'epsc');
    exportgraphics(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.eps')),'BackgroundColor','none','ContentType','vector');
    saveas(gcf,fullfile(mcmcfile, 'figures',strcat('KpolyHistograms_BNI1_NTD',num2str(NTDdist),'.png')),'png');

    % bni1mat=zeros(length(params),Nmax-Nmin+1);
    % 
    % for p=1:length(params)
    %     formin1.opts.k_cap=params(p,1);
    %     formin1.opts.k_del=params(p,2);
    %     formin1.opts.r_cap=params(p,3);
    %     formin1.opts.r_cap_exp=params(p,4);
    %     [ratios,NTD_dists,fh1ength]=NTD_predictions(formin1,Nmax,Nmin);
    %     bni1mat(p,:)=ratios;
    % end
    % 
    % figure;
    % h=histogram2(reshape(bni1mat',1,[]),repmat(NTD_dists,1,length(params)),[30 length(NTD_dists)],'DisplayStyle','tile','ShowEmptyBins','on');
    % a=colorbar;
    % a.Label.String = 'frequency';
    % xlabel("log_{2}(k_{poly} N terminal dimerized/k_{poly} double)")
    % ylabel('Distance from most NT PRM to NTD')
    % 
    % saveas(gcf,fullfile(mcmcresultsloc, 'figures','BNI1kpolyNTDsweep.png'),'png');
    % exportgraphics(gca,fullfile(mcmcresultsloc, 'figures','BNI1kpolyNTDsweep.eps'),'BackgroundColor','none','ContentType','vector');


    function [ratios,NTD_dists,fh1ength]=NTD_predictions(formin,Nmax,Nmin)
    % NTD_predictions sweeps through all possible NTD locations up to fh1 length=Nmax for
    % the input formin and returns kpoly values, corresponding NTD distances
    % and FH1 lengths
    ratios=NaN(1,Nmax-Nmin+1);
    NTD_dists=NaN(1,Nmax-Nmin+1);
    fh1ength=NaN(1,Nmax-Nmin+1);
    for i=0:(Nmax-1)
        if formin.length<Nmax+1
            kpoly=formin.kpoly;
            NTD_dist=formin.PRMList(1,formin.PRMCount).dist_NT;
            ratioval=log2(kpoly.ratio);
            fh1ength(i+1)=formin.length;
            NTD_dists(i+1)=NTD_dist;
            if ratioval==-Inf
                ratioval=NaN;
            end
            ratios(i+1)=ratioval;
            formin.add_length(1)
        else
            formin.add_length(-i)
            return
        end
    end
end
end