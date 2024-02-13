% GENERATE_SWEEP_DATA ...
    % 
    % ...
    %
    % See also .

    % rows:
    % 1) settings
    %


    table_num=233;
    tab=addformin(opts,1,1,1,1);
    tab = table('Size',[3000 length(tab.Properties.VariableNames)],'VariableNames',tab.Properties.VariableNames,'VariableTypes',varfun(@class,tab,'OutputFormat','cell'));
    for FH1_length=100:20:600
        for PRM_size=[1,3,4,5,6,7,9,10,15,20]
            for PRM_loc=[ceil(PRM_size/2),ceil(FH1_length*0.25),ceil(FH1_length/2),ceil(FH1_length*0.75),FH1_length]
                index=0;
                for k=1:1000
                    r=(20--4)*rand(1,5)+-4;
                    opts.k_cap=10^(r(1));
                    opts.k_del=10^(r(2));
                    opts.r_cap=10^(r(3));
                    opts.r_del=10^(r(4));
                    opts.k_rel=10^(r(5));
                    opts.kpoly_type="4st";
                    opts.set_equation(1);

                    index=index+1;
                    tab(index,:)=addformin(opts,1,FH1_length,PRM_loc,PRM_size);

                    opts.set_equation(0,"krel",{"size","negexp"})
                    index=index+1;
                    tab(index,:)=addformin(opts,1,FH1_length,PRM_loc,PRM_size);

                    opts.kpoly_type="3st";
                    index=index+1;
                    tab(index,:)=addformin(opts,1,FH1_length,PRM_loc,PRM_size);
                end
                table_num=table_num+1;
                save(strcat("table_",num2str(table_num),".mat"),"tab",'-mat')
                %savetable(tab,table_num)
                disp(table_num)
            end
        end
    end


function tab=addformin(options,gating,FH1_length,PRM_loc,PRM_size)
    formin1=Formin("formin1",options,c_PA=1,gating=gating,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
    kpoly=formin1.kpoly;
    kcap=formin1.kcap;
    kdel=formin1.kdel;
    rcap=formin1.rcap;
    rdel=formin1.rdel;
    krel=formin1.krel;
    
    tab=table;
    tab.FH1_length=FH1_length;
    tab.PRM_loc=PRM_loc;
    tab.PRM_size=PRM_size;
    tab.gating=gating;
    tab.kpoly_type=formin1.opts.kpoly_type;

    params={"kcap_eq","kdel_eq","rcap_eq","rdel_eq","krel_eq"};
    for i=1:length(params)
        if isfield(formin1.opts.equationstext,params{i})
            tab.(params{i})=formin1.opts.equationstext.(params{i});
        else
            tab.(params{i})="n/a";
        end
    end

    params={"k_cap","k_del","r_cap","r_del","k_rel"};
    for i=1:length(params)
        tab.(params{i})=formin1.opts.(params{i});
    end

    rates={kpoly,kcap,kdel,rcap,rdel,krel};
    params={"kpoly","kcap","kdel","rcap","rdel","krel"};
    params_ratio={"kpoly_ratio","kcap_ratio","kdel_ratio","rcap_ratio","rdel_ratio","krel_ratio"};
    params_double={"kpoly_double","kcap_double","kdel_double","rcap_double","rdel_double","krel_double"};
    params_dimer={"kpoly_dimer","kcap_dimer","kdel_dimer","rcap_dimer","rdel_dimer","krel_dimer"};
    for i=1:length(rates)
        if class(rates{i})=="FilType"
            tab.(params_ratio{i})=rates{i}.ratio;
            tab.(params_double{i})=rates{i}.double;
            tab.(params_dimer{i})=rates{i}.dimer;
        else
            tab.(params{i})=rates{i};
        end
    end

    stats={"Prvec0","POcclude","POcclude_Base"};
    stats_ratio={"Prvec0_ratio","POcclude_ratio","POcclude_Base_ratio"};
    stats_doublea={"Prvec0_doublea","POcclude_doublea","POcclude_Base_doublea"};
    stats_doubleb={"Prvec0_doubleb","POcclude_doubleb","POcclude_Base_doubleb"};
    stats_dimera={"Prvec0_dimera","POcclude_dimera","POcclude_Base_dimera"};
    stats_dimerb={"Prvec0_dimerb","POcclude_dimerb","POcclude_Base_dimerb"};
    for i=1:length(stats)
        stat=formin1.PRMList(1).(stats{i});
        tab.(stats_ratio{i})=stat.ratio;
        tab.(stats_doublea{i})=stat.double.a;
        tab.(stats_doubleb{i})=stat.double.b;
        tab.(stats_dimera{i})=stat.dimer.a;
        tab.(stats_dimerb{i})=stat.dimer.b;
    end

end

function savetable(table,num)
    table_temp=table;
    save(strcat("table_",num2str(num),".mat"),"table_temp",'-mat')
end
