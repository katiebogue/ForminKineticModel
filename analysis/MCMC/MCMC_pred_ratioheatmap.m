function MCMC_pred_ratioheatmap(E1,E2, E3,titles,R1, R2, R3)
%MCMC_PRED_RATIOHEATMAP 
%
%   Inputs:
%       E1 : matrix of ratios
%       E2 : matrix of ratios
%       E2 : matrix of ratios
%       titles: 1x3 list of strings
%       R1 : 1x2 matrix where R1(1) is the minimum bounds of the error bars
%       for E1 and R1(2) is the maximum bounds (0,0 if no bounds)
%       R2 : 1x2 matrix where R2(1) is the minimum bounds of the error bars
%       for E2 and R2(2) is the maximum bounds (0,0 if no bounds)
%       R3 : 1x2 matrix where R3(1) is the minimum bounds of the error bars
%       for E3 and R3(2) is the maximum bounds (0,0 if no bounds)
%
%   
%
%   Runs trueparamconvert.
%
% See also TRUEPARAMCONVERT, MCMC_PRED.
    NBINS = 30;
    figure('units','centimeters','position',[5,5,35,15],'Name','Particle clouds');hold on;
    t=tiledlayout(1,3,'TileSpacing','tight','Padding','none');
    t.Title.String = 'Kpoly ratios (dimer/double)';

    nexttile(1); hold on;
    hist3([E1 E2],[NBINS NBINS],'CdataMode','auto','LineStyle','none')
    c=colormap('cool');
    colormap([1 1 1; c])
    view(2)
    xlabel(titles(1))
    ylabel(titles(2))
    if any(R2)
        yline(R2(1),'--r','LineWidth',2)
        yline(R2(2),'--r','LineWidth',2)
    end
    if any(R1)
        xline(R1(1),'--r','LineWidth',2)
        xline(R1(2),'--r','LineWidth',2)
    end

    nexttile(2); hold on;
    hist3([E1 E3],[NBINS NBINS],'CdataMode','auto','LineStyle','none')
    c=colormap('cool');
    colormap([1 1 1; c])
    view(2)
    xlabel(titles(1))
    ylabel(titles(3))
    if any(R3)
        yline(R3(1),'--r','LineWidth',2)
        yline(R3(2),'--r','LineWidth',2)
    end
    if any(R1)
        xline(R1(1),'--r','LineWidth',2)
        xline(R1(2),'--r','LineWidth',2)
    end

    nexttile(3); hold on;
    hist3([E3 E2],[NBINS NBINS],'CdataMode','auto','LineStyle','none')
    c=colormap('cool');
    colormap([1 1 1; c])
    view(2)
    xlabel(titles(3))
    ylabel(titles(2))
    if any(R2)
        yline(R2(1),'--r','LineWidth',2)
        yline(R2(2),'--r','LineWidth',2)
    end
    if any(R3)
        xline(R3(1),'--r','LineWidth',2)
        xline(R3(2),'--r','LineWidth',2)
    end

    a=colorbar;
    a.Label.String = 'Frequency';
	fontsize(15,"points")
end