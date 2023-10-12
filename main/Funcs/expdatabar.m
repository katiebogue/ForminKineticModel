function fig = expdatabar(datastruct)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    arguments
        datastruct struct
    end
    
    formins=[datastruct.formin];
    kpolys=[formins.kpoly];
    types=[datastruct.type];

    simdata=zeros(1,length(datastruct));
    
    for i=1:length(datastruct)
        type=types(i);
        simdata(i)=kpolys(i).(type);
    end

    fig=figure().empty;
    makeplot("all");
    makeplot("single");
    makeplot("double");
    makeplot("dimer");
    makeplot("ratio");

    function makeplot(type)
        if type=="all"
            typesimdata=simdata;
            typevalues=datastruct;
        else
            typesimdata=simdata(types==type);
            typevalues=datastruct(types==type);
        end

        if length(typesimdata)==0
            return
        end
        
        fig(1+end)=figure;
        
        if length(typesimdata)==1
            kpolybar=bar(1,[[typevalues.value]',typesimdata']);
        else
            kpolybar=bar([[typevalues.value]',typesimdata']);
        end

        typeformins=[typevalues.formin];
    
        set(gca,'xticklabel',[typeformins.name]);
        hold on
        
        erlocs=(0:(length(typesimdata)-1))+0.85;
        er = errorbar(erlocs,[typevalues.value]', [typevalues.errbot],[typevalues.errtop]);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        
        set(kpolybar(1), 'FaceColor','b');
        set(kpolybar(2), 'FaceColor','r');
        legend( 'Experimental', 'model');
        xlabel('Formins');
        ylabel(strcat("Kpoly ",type));
    end
end