function fig=varltplot(obj,xval,stat,skip)
    arguments
        obj Lookuptable
        xval string {mustBeMember(xval,{'length','NTdist','CTdist'})}
        stat string %need to figure out validation
        skip double =1 % will plots points from 1:skip:end
    end
    obj.holdratio=true;

    hold on
    fig=figure;

    mat_dimer=obj.stattable(stat,"dimer");
    mat_double=obj.stattable(stat,"double");
    %mat_ratio=obj.stattable(stat,"ratio");
    mat_dimer_diff=mat_dimer.a;
    mat_dimer_diff(:,3)=abs(mat_dimer.a(:,3)-mat_dimer.b(:,3));
    mat_double_diff=mat_double.a;
    mat_double_diff(:,3)=abs(mat_double.a(:,3)-mat_double.b(:,3));
    mat_diffa=mat_double.a;
    mat_diffa(:,3)=abs(mat_double.a(:,3)-mat_dimer.a(:,3));
    mat_diffb=mat_double.b;
    mat_diffb(:,3)=abs(mat_double.b(:,3)-mat_dimer.b(:,3));

    mat_diffa_dimer=mat_double.b;
    mat_diffa_dimer(:,3)=mat_diffa(:,3)-mat_dimer_diff(:,3);
    indx=mat_diffa_dimer(:,3)>0;
    mat_diffa_dimer(indx,3)=0;
    if class(mat_dimer)=="Filament"
        %filscatter(mat_ratio.a,"ratio a",'#29ABE2',"none")
        %filscatter(mat_ratio.b,"ratio b",'#F15A22',"none")
        filscatter(mat_diffb,"double-dimer b",'#F15A22',"none")
        filscatter(mat_diffa,"double-dimer a",'#29ABE2',"none")
        filscatter(mat_dimer_diff,"dimer diff",'red',"none")
        filscatter(mat_double_diff,"double diff",'blue',"none")

        filscatter(mat_diffa_dimer,"double-dimer a -dimer diff",'black',"none")
        legend
    else
        error("stat must be of type filament")
    end

    function filscatter(mat,label,color,scale)
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
        % lne=yline(ratioline,'LineWidth',2);
        % lne.Annotation.LegendInformation.IconDisplayStyle='off';
        s.MarkerFaceColor=color;

        xlabel(xlab)
        ylabel(ylab)
    end

    hold off
    obj.holdratio=false;
end