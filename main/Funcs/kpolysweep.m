function fig=kpolysweep(options,sweep,scale,save)
% KPOLYSWEEP  ...
    %
    %   KPOLYSWEEP() ...
    %
    %   Inputs:
    %       formin  : 
    %
    %
    %   See also .

arguments
    options Options
    sweep string
    scale string {mustBeMember(scale,{'none','log2','log10','ln'})}="log2"
    save logical=false
end
fig=figure;
hold on

if scale=="none"
    yscale=@(x) x;
    ylab="k_{poly} N terminal dimerized/k_{poly} double";
elseif scale=="log2"
    ylab="log_{2}(k_{poly} N terminal dimerized/k_{poly} double)";
    yscale=@(x) log2(x);
elseif scale=="log10"
    ylab="log_{10}(k_{poly} N terminal dimerized/k_{poly} double)";
    yscale=@(x) log10(x);
elseif scale=="ln"
    ylab="ln(k_{poly} N terminal dimerized/k_{poly} double)";
    yscale=@(x) log(x);
end

colors=["blue","red","black"];

if sweep=="PRM size"
    FH1_lengths=[50,150,350,600];
    scatn=0;
    for j=1:length(FH1_lengths)
        FH1_length=FH1_lengths(j);
        for k=0.1:0.4:0.9
            if k==.1
                color1='blue';
            elseif k==0.5
                color1='red';
            else
                color1='black';
            end

            if FH1_length<60 && k>0.5
                continue
            end
            maxval=min(50,FH1_length);
            xvals=[1:2:10,10:10:maxval];
            kpolyratios=zeros(length(xvals),1);
            PRM_loc=ceil(FH1_length*k);
            for i=1:length(xvals)
                formin1=Formin("formin1",options,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=xvals(i));
                kpolyratios(i,1)=formin1.kpoly.ratio;
            end
            scatn=scatn+1;
            %scat=scatter(xvals',yscale(kpolyratios),'filled',options.shapes(j),'MarkerFaceColor',options.colors(scatn));
            scat=scatter(xvals',yscale(kpolyratios),'filled',options.shapes(j),'MarkerFaceColor',color1);
            %plt=plot(xvals',yscale(kpolyratios),'Color',options.colors(scatn));
            plt=plot(xvals',yscale(kpolyratios),'Color',color1);
            plt.Annotation.LegendInformation.IconDisplayStyle='off';
            scat.DisplayName=sprintf("FH1 length=%d; FH2 dist=%d",FH1_length,PRM_loc);
        end
    end
    xlab="PRM size";
    xlabel(xlab)
elseif sweep=="FH1 length"
    PRM_sizes=[4,7,10];
    scatn=0;
    for j=1:length(PRM_sizes)
        PRM_size=PRM_sizes(j);
        PRM_locs=[10,50,100,200];
        for k=1:length(PRM_locs)
            PRM_loc=PRM_locs(k);
            xvals=[PRM_loc:10:400];
            kpolyratios=zeros(length(xvals),1);
            for i=1:length(xvals)
                formin1=Formin("formin1",options,c_PA=1,gating=1,length=xvals(i),PRMloc=PRM_loc,PRMsize=PRM_size);
                kpolyratios(i,1)=formin1.kpoly.ratio;
            end
            scatn=scatn+1;
            scat=scatter(xvals',yscale(kpolyratios),'filled',options.shapes(k),'MarkerFaceColor',colors(j));
            plt=plot(xvals',yscale(kpolyratios),'Color',colors(j));
            % scat=scatter(xvals',yscale(kpolyratios),'filled',options.shapes(k),'MarkerFaceColor',options.colors(scatn));
            % plt=plot(xvals',yscale(kpolyratios),'Color',options.colors(scatn));
            plt.Annotation.LegendInformation.IconDisplayStyle='off';
            scat.DisplayName=sprintf("PRM size=%d; FH2 dist=%d",PRM_size,PRM_loc);
        end
    end
    xlab="FH1 length";
    xlabel(xlab)
elseif sweep=="FH2 dist"
    FH1_lengths=[50,150,300,400];
    scatn=0;
    for j=1:length(FH1_lengths)
        FH1_length=FH1_lengths(j);
        PRM_sizes=[4,7,10];
        for k=1:length(PRM_sizes)
            xvals=1:5:FH1_length;
            PRM_size=PRM_sizes(k);
            kpolyratios=zeros(length(xvals),1);
            for i=1:length(xvals)
                PRM_loc=xvals(i);
                formin1=Formin("formin1",options,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
                kpolyratios(i,1)=formin1.kpoly.ratio;
            end
            scatn=scatn+1;
            scat=scatter(xvals',yscale(kpolyratios),'filled',options.shapes(j),'MarkerFaceColor',colors(k));
            plt=plot(xvals',yscale(kpolyratios),'Color',colors(k));
            plt.Annotation.LegendInformation.IconDisplayStyle='off';
            scat.DisplayName=sprintf("FH1 length=%d; PRM size=%d",FH1_length,PRM_size);
        end
    end
    xlab="Distance from PRM to FH2";
    xlabel(xlab)
elseif sweep=="NT dist"
    FH1_lengths=[50,150,300,400];
    scatn=0;
    for j=1:length(FH1_lengths)
        FH1_length=FH1_lengths(j);
        PRM_sizes=[4,7,10];
        for k=1:length(PRM_sizes)
            xvals=FH1_length:-5:1;
            PRM_size=PRM_sizes(k);
            kpolyratios=zeros(length(xvals),1);
            for i=1:length(xvals)
                PRM_loc=xvals(i);
                formin1=Formin("formin1",options,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
                kpolyratios(i,1)=formin1.kpoly.ratio;
            end
            scatn=scatn+1;
            scat=scatter([FH1_length-xvals]',yscale(kpolyratios),'filled',options.shapes(j),'MarkerFaceColor',colors(k));
            plt=plot([FH1_length-xvals]',yscale(kpolyratios),'Color',colors(k));
            plt.Annotation.LegendInformation.IconDisplayStyle='off';
            scat.DisplayName=sprintf("FH1 length=%d; PRM size=%d",FH1_length,PRM_size);
        end
    end
    xlab="Distance from PRM to N-term";
    xlabel(xlab)
else
    if ismember(sweep,fieldnames(options))
        ogopt=options.(sweep);
        FH1_lengths=[50,150,500];
        PRM_sizes=[4,10];
        scatn=0;
        for j=1:length(FH1_lengths)
            FH1_length=FH1_lengths(j);
            for k=[0.25,0.75]
                if FH1_length<60 && k>0.5
                    %continue
                end
                PRM_loc=ceil(FH1_length*k);
                for i=1:length(PRM_sizes)
                    PRM_size=PRM_sizes(i);
                    formin1=Formin("formin1",options,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
                    if sweep=="k_cap"
                        yvals=10.^([0:0.2:10]);
                    elseif sweep=="k_del"
                        yvals=10.^([-8:0.2:6]);
                    elseif sweep=="r_cap"
                        yvals=10.^([-2:0.2:12]);
                    elseif sweep=="r_del"
                        yvals=10.^([0:0.2:10]);
                    elseif sweep=="k_rel"
                        yvals=10.^([-2:0.2:14]);
                    elseif sweep=="r_cap_exp"
                        yvals=10.^([-2:0.1:2]);
                    else
                        yvals=10.^([-2:0.2:10]);
                    end
                    kpolyratios=zeros(length(yvals),1);
                    for n=1:length(yvals)
                        options.(sweep)=yvals(n);
                        kpolyratios(n,1)=formin1.kpoly.ratio;
                    end
                    scatn=scatn+1;
                    if k==0.75
                        scat=scatter(log10(yvals),yscale(kpolyratios),options.shapes(j),'MarkerEdgeColor',colors(i));
                    else
                        scat=scatter(log10(yvals),yscale(kpolyratios),15,'filled',options.shapes(j),'MarkerFaceColor',colors(i));
                    end
                    plt=plot(log10(yvals),yscale(kpolyratios),'Color',colors(i));
                    % scat=scatter(log10(yvals),yscale(kpolyratios),'filled',options.shapes(j),'MarkerFaceColor',options.colors(scatn));
                    % plt=plot(log10(yvals),yscale(kpolyratios),'Color',options.colors(scatn));
                    plt.Annotation.LegendInformation.IconDisplayStyle='off';
                    scat.DisplayName=sprintf("FH1 length=%d; FH2 dist=%d; PRM size=%d",FH1_length,PRM_loc,PRM_size);
                end
            end
        end
        options.(sweep)=ogopt;
        strchar=char(sweep);
        if strchar(1:2)=="k_" || strchar(1:2)=="r_"
            sweep=strcat(strchar(1:2),"{",strchar(3:end),"}");
        end
        xlab=strcat("log_{10}(",sweep,")");
        xlabel(xlab)
    else
        error("please enter a valid sweep method")
    end
end

ylabel(ylab)
title(strcat("Change in elongation rate vs. ",xlab))
legend('Location','southoutside','NumColumns',2,'FontSize',8)
lne=yline(yscale(1),'LineWidth',2);
lne.Annotation.LegendInformation.IconDisplayStyle='off';

if save
    figuresave(gcf,options,append("sweep_",gca().Title.String,'.fig'),true);
end

end