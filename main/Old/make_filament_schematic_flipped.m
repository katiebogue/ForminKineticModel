% can use in place of make_filiment_schematic

% outputs line and circle graph displaying binding site location (circles)
% across the entire length of an fh1 filament (line) 
% as in make_filament_schematic

% graphs filament schematics with x and y axes switched;
% i.e. horizantal filament schematic

% uses variables: pp_length_vec, pp_index_vec, iSite_tot, fh1_length,
% fh1_name, opt1/2/3/4

%sets formatting options based on opt4
if opt4 == 1
    if LOOP == 1
        figure()
        subplot(3,4,1)
    else
        if rem(LOOP,3) == 1
            figure()
            subplot(3,4,1)
        elseif rem(LOOP,3) == 2
            subplot(3,4,5)
        elseif rem(LOOP,3) == 0
            subplot(3,4,9)
        end
    end
    
else
    figure()
    subplot(1,4,1)
end

hold on
x = [1,1];
y = [1,fh1_length];

plot(y,x,'k','LineWidth',2)

for i = 1:iSite_tot
    location = pp_index_vec(i);
    width = pp_length_vec(i);
    scatter(location,1,12*width,'filled')
end

set(gca,'YTick',[])
if opt4 == 0
    xlabel('C-Terminus     --     Amino Acid Index     --     N-Terminus')
else
    xlabel('C -- N-Terminus')
end
ylabel('FH1')

if opt == 1
    xlim([1,fh1_length])
elseif opt ==2
    if fh1_length <= 300
        xlim([1,300])
    else
        xlim([1,fh1_length])
    end  
else
    xlim([1,515])
end

hold off

if opt2 == 1
    if opt4 == 1
        title(fh1_name);
    else
        sgtitle(fh1_name);
    end
else
    if opt4 == 1
        title(fh1_info);
    else
        sgtitle(fh1_info);
    end
end