function fig=ltplot(obj,xval,stat,skip,NameValueArgs)
    arguments
        obj Lookuptable
        xval string {mustBeMember(xval,{'length','NTdist','CTdist'})}
        stat string %need to figure out validation
        skip double =1 % will plots points from 1:skip:end
        NameValueArgs.type string {mustBeMember(NameValueArgs.type,{'single','dimer','double','ratio'})}
        NameValueArgs.ratioscale string {mustBeMember(NameValueArgs.ratioscale,{'none','log2','log10','ln'})} % scale for ratio plot
        NameValueArgs.scale string {mustBeMember(NameValueArgs.scale,{'none','log2','log10','ln'})} % will set scale for all plots (including ratio if no ratioscale is set)
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
        fig(1)=figure;
        if NameValueArgs.type=="ratio"
            typscatter(NameValueArgs.type,ratioscale)
        else
            typscatter(NameValueArgs.type,scale)
        end
        legend
    else
        fig(1)=figure;
        typscatter("single",scale)
        typscatter("double",scale)
        typscatter("dimer",scale)
        legend
        fig(2)=figure;
        typscatter("ratio",ratioscale)
        legend
    end

    hold off
    obj.holdratio=false;

    function typscatter(type,scale)
        mat=obj.stattable(stat,type);
        if class(mat)=="Filament"
            filscatter(mat.a,strcat(type," a"),'#29ABE2')
            filscatter(mat.b,strcat(type," b"),'#F15A22')
        else
            filscatter(mat,type,'#F15A22')
        end

        function filscatter(mat,label,color)
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
            
            
            s=scatter(x(1:skip:end),y(1:skip:end), 'filled','DisplayName',label);
            if type=="ratio"
                lne=yline(ratioline,'LineWidth',2);
                lne.Annotation.LegendInformation.IconDisplayStyle='off';
            end
            s.MarkerFaceColor=color;
    
            xlabel(xlab)
            ylabel(ylab)
        end
    end


end