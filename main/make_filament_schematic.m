% outputs line and circle graph displaying binding site location (circles)
% across the entire length of an fh1 filament (line)

% uses variables: pp_length_vec, pp_index_vec, iSite_tot, fh1_length,
% fh1_name, opt1/2/3/4

%sets formatting options based on opt4
if opt4 == 1
    if LOOP == 1 %first formin
        figure()
        subplot(3,4,1) %3x4 grid w/ axis on 1st cell
    else
        if rem(LOOP,3) == 1  %divisible by 3 w/ remainder 1
            figure()
            subplot(3,4,1) %3x4 grid w/ axis on 1st cell
        elseif rem(LOOP,3) == 2  %divisible by 3 w/ remainder 2
            subplot(3,4,5) %3x4 grid w/ axis on 5th cell
        elseif rem(LOOP,3) == 0 %if divisible by 3
            subplot(3,4,9) %3x4 grid w/ axis on 9th cell
        end
    end
    
else
    figure()
    subplot(1,4,1) %1x4 grid w/ axis on 1st cell
end

hold on
x = [1,1];
y = [1,fh1_length];

%plots line to represent fh1 chain
plot(x,y,'k','LineWidth',2)

% plots circles to represent binding sites
% width of circle depends on number of polyprolines at the binding site
for i = 1:iSite_tot
    location = pp_index_vec(i);
    width = pp_length_vec(i);
    scatter(1,location,12*width,'filled') %adds dots with width proportional to # of PP
end

set(gca,'XTick',[])
if opt4 == 0
    ylabel('C-Terminus     --     Amino Acid Index     --     N-Terminus')
else
    ylabel('C -- N-Terminus')
end
xlabel('FH1')

% formatting options based on opt1 and opt2
if opt1 == 1
    ylim([1,fh1_length])
elseif opt1 ==2
    if fh1_length <= 300
        ylim([1,300])
    else
        ylim([1,fh1_length])
    end  
else
    ylim([1,515])
end

hold off

if opt2 == 1
    if opt4 == 1
        title(fh1_name); %just name
    else
        sgtitle(fh1_name);
    end
else
    if opt4 == 1
        title(fh1_info); %all info in fh1_info
    else
        sgtitle(fh1_info);
    end
end


