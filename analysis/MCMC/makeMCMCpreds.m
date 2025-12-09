function makeMCMCpreds(mcmcresultsloc)
    load("forminexperimentobjs.mat","Experiment_BNI1")
    load("forminexperimentobjs.mat","Experiment_FHOD")
    load("forminexperimentobjs.mat","Experiment_Capu")

    kpoly_ratios_FHOD=MCMC_pred(Experiment_FHOD, fullfile(mcmcresultsloc, "mcmc_results.mat"),0, 0.6042, 0.6996);
    sgtitle("FHOD")
    %saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_FHOD.eps'),'epsc');
    saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_FHOD.png'),'png');
    
    kpoly_ratios_BNI1=MCMC_pred(Experiment_BNI1, fullfile(mcmcresultsloc, "mcmc_results.mat"), 0, 0, 0);
    sgtitle("BNI1")
    %saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_BNI1.eps'),'epsc');
    saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_BNI1.png'),'png');

    kpoly_ratios_CAPU=MCMC_pred(Experiment_Capu, fullfile(mcmcresultsloc, "mcmc_results.mat"),0, 0.9523, 1.0477);
    sgtitle("CAPU")
    %saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_CAPU.eps'),'epsc');
    saveas(gcf,fullfile(mcmcresultsloc, 'figures','KpolyHistograms_CAPU.png'),'png');

    MCMC_pred_ratioheatmap(kpoly_ratios_BNI1,kpoly_ratios_CAPU,kpoly_ratios_FHOD,["BNI1","CAPU","FHOD"], [0 0],[0.9523 1.0477],[0.6042 0.6996])
    saveas(gcf,fullfile(mcmcresultsloc, 'figures','ratioheatmapsBNI1_CAPU_FHOD.png'),'png');
end