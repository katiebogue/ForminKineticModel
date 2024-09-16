function fig=ltplot(obj,xval,stat,skip,NameValueArgs)
    arguments
        obj Lookuptable
        xval string {mustBeMember(xval,{'length','NTdist','CTdist'})}
        stat string %need to figure out validation
        skip double =1 % will plots points from 1:skip:end
        NameValueArgs.type string {mustBeMember(NameValueArgs.type,{'single','dimer','double','ratio'})}
        NameValueArgs.ratioscale string {mustBeMember(NameValueArgs.ratioscale,{'none','log2','log10','ln'})} % scale for ratio plot
        NameValueArgs.scale string {mustBeMember(NameValueArgs.scale,{'none','log2','log10','ln'})} % will set scale for all plots (including ratio if no ratioscale is set)
        NameValueArgs.ax1 %axes to put the plot on
    end
    obj.holdratio=true;

    if isfield(NameValueArgs,"scale")
        scale=NameValueArgs.scale;
        if isfield(NameValueArgs,"ratioscale")
            ratioscale=NameValueArgs.ratioscale;
        else
            ratioscale=scale;
        end
    else
        scale="none";
        if isfield(NameValueArgs,"ratioscale")
            ratioscale=NameValueArgs.ratioscale;
        else
            ratioscale="log2";
        end
    end

    if isfield(NameValueArgs,"type")
        if isfield(NameValueArgs,"ax1")
            ax1=NameValueArgs.ax1;
        else
            fig(1)=figure;
            ax1=gca;
        end
        if NameValueArgs.type=="ratio"
            typscatter(NameValueArgs.type,ratioscale,ax1)
        else
            typscatter(NameValueArgs.type,scale,ax1)
        end
        legend
    else
        fig(1)=figure;
        typscatter("single",scale,gca)
        typscatter("double",scale,gca)
        typscatter("dimer",scale,gca)
        legend
        fig(2)=figure;
        typscatter("ratio",ratioscale,gca)
        legend
    end

    set(gca,'fontname','Arial')
    hold off
    obj.holdratio=false;

    function typscatter(type,scale,ax1)
        mat=obj.stattable(stat,type);
        if class(mat)=="Filament"
            filscatter(mat.a,strcat(type," a"),'#29ABE2',ax1)
            filscatter(mat.b,strcat(type," b"),'#F15A22',ax1)
        else
            filscatter(mat,type,'#F15A22',ax1)
        end

        function filscatter(mat,label,color,ax1)
            hold on
            if xval=="length"
                x=mat(:,1);
                xlab="FH1 length";
            elseif xval=="NTdist"
                x=mat(:,1)-mat(:,2);
                xlab="Distance from PRM to N-terminus";
            elseif xval=="CTdist"
                x=mat(:,2);
                xlab="Distance from PRM to FH2";
            end
            
            y=mat(:,3);
            ylab=strcat(scale,"(",stat,")");

            if scale=="none"
                ylab=stat;
                ratioline=1;
            elseif scale=="log2"
                y=log2(y);
                ratioline=0;
            elseif scale=="log10"
                y=log10(y);
                ratioline=0;
            elseif scale=="ln"
                y=log(y);
                ratioline=0;
            end
            
            
            s=scatter(ax1,x(1:skip:end),y(1:skip:end), 'filled','DisplayName',label);
            if type=="ratio"
                lne=yline(ax1,ratioline,'LineWidth',1,'Color',[0 0 0]);
                lne.Annotation.LegendInformation.IconDisplayStyle='off';
            end
            s.MarkerFaceColor=color;
    
            xlabel(ax1,xlab)
            ylabel(ax1,ylab)
        end
    end


end