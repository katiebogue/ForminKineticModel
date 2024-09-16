function fig = filamentSchematic(formin)
% FILAMENTSCHEMATIC  creates schematic of the formin filament.
    %   fig = FILAMENTSCHEMATIC(formin) creates schematic of formin. 
    %
    %   Schematic is 2 lines representing the entire length of an fh1
    %   filament. 1 line has a colored dot for the center of each binding 
    %   site; 1 line has blue dots for each proline, with binding sites 
    %   highlighted in yellow. 
    %
    %   Inputs:
    %       formin : Formin object
    %   
    %   See also FORMIN, PRM, EXPERIMENT, OPTIONS.

    arguments
        formin Formin
    end
    % outputs line and circle graph displaying binding site location (circles)
    % across the entire length of an fh1 filament (line)
    
    
    FH1_length=formin.length;
    PRMCount=formin.PRMCount;
    PRMLoc=[formin.PRMList.dist_FH2];
    PRMLoc_start=[formin.PRMList.dist_FH2_start];
    PRMLen=[formin.PRMList.size];
    PLoc=formin.Ploc;
    colors=formin.opts.colors;

    fig=figure;
    hold on
    x = [1,1];
    z = [1.5,1.5];
    y = [1,FH1_length];
    
    
    %plots line to represent fh1 chain
    plot(x,y,'k','LineWidth',2)
    plot(z,y,'k','LineWidth',2)
    
    
    % plots circles to represent binding sites
    % width of circle depends on number of polyprolines at the binding site
    for i = 1:PRMCount
        location = PRMLoc(i);
        width = PRMLen(i);
        scatter(1,location,12*width,'filled','MarkerFaceColor', colors(i)) %adds dots with width proportional to # of PP
    end
    for i = 1:PRMCount
        start = PRMLoc_start(i)-0.6; 
        width = PRMLen(i);
        eend = start+width+0.4;
        ly = [start,eend];
        plot(z,ly,'color','#EDB120','LineWidth',8) %adds highlight of PRMS
    end
    for i = 1:length(PLoc)
        location = PLoc(i);
        width = 2;
        scatter(1.5,location,12*width,'filled','MarkerFaceColor','#0072BD','MarkerEdgeColor','k') %adds dots for each P
    end
    
    xlim([0.5,2])
    
    set(gca,'XTick',[])
    ylabel('C-Terminus     --     Amino Acid Index     --     N-Terminus')
    xlabel('FH1')

    fig.Position = [455 354 150 506];
    
end