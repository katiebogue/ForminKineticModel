function fig=ltplot(obj,xval,stat,skip,NameValueArgs)
% LTPLOT create scatterplots of the specified polymer stat vs. an FH1 property
% 
% fig = LOOKUPTABLE.LTPLOT(xval,stat,skip,NameValueArgs)
%
%   Inputs: 
%       xval        : (string) FH1 property to scatter points by (options:
%                   'length','NTdist','CTdist')
%       stat        : (string) polymer stat to plot, must be a lookuptable
%                   property
%       skip        : (double) plot points from 1:skip:end to improve plot
%                   readability (default is 1-- aka no skip)
%       type        : (string) "dimer," "double," "single," or "ratio" (NameValueArgs)
%       FH1         : (double) size of the FH1 to fix; only works properly if xval='NTdist' or 'CTdist' (use -1 to plot all;
%                    default is -1) (NameValueArgs)
%       PRM         : (double) size of the FH1 to fix; only works properly if xval='length' (use -1 to plot all;
%                    default is -1) (NameValueArgs)
%       ratioscale  : (string) how to scale ratio y axis (none, log2,
%                   log10,ln), if not specified, used the same scale as
%                   specifed by scale (NameValueArgs) 
%       scale       : (string) how to scale stat property y axis (none, log2,
%                   log10,ln, amino((1.0e33*var/27*6.022e23)), if not specified, no scaling is applied (NameValueArgs)
%       ax1         : (axes) axes to put the plot on (NameValueArgs)
% 
%   Output is an array of figures holding the scatterplots.
% 
% See also LOOKUPTABLE, LOOKUPTABLE/STATTABLE.
    arguments
        obj Lookuptable
        xval string {mustBeMember(xval,{'length','NTdist','CTdist'})}
        stat string 
        skip double =1 % will plots points from 1:skip:end
        NameValueArgs.type string {mustBeMember(NameValueArgs.type,{'single','dimer','double','ratio'})}
        NameValueArgs.FH1 double
        NameValueArgs.PRM double
        NameValueArgs.ratioscale string {mustBeMember(NameValueArgs.ratioscale,{'none','log2','log10','ln'})} % scale for ratio plot
        NameValueArgs.scale string {mustBeMember(NameValueArgs.scale,{'none','log2','log10','ln','amino'})} % will set scale for all plots (including ratio if no ratioscale is set)
        NameValueArgs.ax1 %axes to put the plot on
    end
    obj.holdratio=true;

    fixFH1=false;
    if isfield(NameValueArgs,"FH1")
        if NameValueArgs.FH1~=-1
            fixFH1=true;
            FH1fixedval=NameValueArgs.FH1;
        end
    end

    fixPRM=false;
    if isfield(NameValueArgs,"PRM")
        if NameValueArgs.PRM~=-1
            fixPRM=true;
            PRMfixedval=NameValueArgs.PRM;
        end
    end

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
        if fixFH1
            mat.a = mat.a(mat.a(:, 1) == FH1fixedval, :); % Filter by fixed FH1 value
            mat.b = mat.b(mat.b(:, 1) == FH1fixedval, :); % Filter by fixed FH1 value
        end
        if fixPRM
            mat.a = mat.a(mat.a(:, 2) == PRMfixedval, :); % Filter by fixed PRM value
            mat.b = mat.b(mat.b(:, 2) == PRMfixedval, :); % Filter by fixed PRM value
        end
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
            elseif scale=="amino"
                y=(1.0e33.*y./(27.*6.022e23));
                ratioline=1;
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